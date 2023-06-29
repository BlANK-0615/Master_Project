import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['NUMEXPR_MAX_THREADS'] = '1'

# import other libraries
import sys
import math
import json
import psutil
import shutil
import pickle
import logging
import argparse
import textwrap
import contextlib
import numpy as np
import pandas as pd
from tqdm import tqdm
import xgboost as xgb
import tensorflow as tf
from utils.ecifs import *
from functools import reduce
from utils import binana, kier
tf.get_logger().setLevel('ERROR')
from utils.dock_functions import *
from functools import partialmethod
from warnings import filterwarnings
from tensorflow.keras.models import load_model
from joblib import Parallel, parallel, delayed, load

#######################################################################
# Global Variables

# filter pandas warnings
filterwarnings('ignore')

# get working directory where scoring function is being deployed
stem_path = os.getcwd()

#######################################################################
# Functions

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):

    """
    Context manager to allow joblib
    progress monitoring
    """

    def tqdm_print_progress(self):
        if self.n_completed_tasks > tqdm_object.n:
            n_completed = self.n_completed_tasks - tqdm_object.n
            tqdm_object.update(n=n_completed)

    original_print_progress = parallel.Parallel.print_progress
    parallel.Parallel.print_progress = tqdm_print_progress

    try:
        yield tqdm_object
    finally:
        parallel.Parallel.print_progress = original_print_progress
        tqdm_object.close()





def prepare_and_dock_inputs(params):

    """
    Function: Prepare smiles inputs for scoring by converting
              to pdbqt files and docking with GWOVina

    Inputs:   Parsed command line parameters

    Outputs:  Updated parameters where the ['ligand'] is docked pdbqt versions
              of supplied smiles strings
              Dictionary of smiles strings with their identifiers
    """

    # load the docking settings
    dock_settings = json.load(open(os.path.join('utils','params','dock_settings.json')))

    # check there has been a binding site location passed to use for docking
    if params.ref_lig is None:
        if params.center is None and params.range is None:

            # if not then complain and exit
            logging.critical("ERROR: No reference ligand or binding site coordinates supplied. Try:\n- ensuring center and range values are entered correctly\n- supplying a reference ligand")
            sys.exit()

        # if center and range have been supplied then parse them into
        # coordinates variable
        else:
            try:
                center_coords = json.loads(params.center)
                center_range = json.loads(params.range)
                coords = (float(center_coords[0]),
                            float(center_coords[1]),
                            float(center_coords[2]),
                            float(center_range[0]),
                            float(center_range[1]),
                            float(center_range[2]))
            
            # if this doesn't work then complain and exit
            except:
                logging.critical("\nERROR: Binding site coordinates for docking are missing or incorrectly defined. \nTry:\n- ensuring center and range values are entered correctly\n- using a reference ligand instead")
                sys.exit()

    # if the user has passed a reference ligand 
    else:
        coords = get_coordinates(params.ref_lig, dock_settings['padding'])

    # make sure all temporary folders exist
    if not os.path.isdir(os.path.join('utils','temp','pdb_files')):
        os.makedirs(os.path.join('utils','temp','pdb_files'))
        os.makedirs(os.path.join('utils','temp','pdbqt_files'))
        os.makedirs(os.path.join('utils','temp','docked_pdbqt_files'))

    # make sure the folders are empty if they exist
    pdbs = get_filepaths(os.path.join('utils','temp','pdb_files',''))
    for pdb in pdbs:
        os.remove(pdb)

    pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))
    for pdbqt in pdbqts:
        os.remove(pdbqt)

    docked_pdbqts = get_filepaths(os.path.join('utils','temp','docked_pdbqt_files',''))
    for docked_pdbqt in docked_pdbqts:
        os.remove(docked_pdbqt)

    # load the smiles from the input file as a dictionary
    # with an identifier as the key and the smile string as the value
    smi_dict = get_smiles(params.ligand)

    logging.info('Generating 3D pdbs from SMILES...')

    # parallelise the conversion of smiles strings to 3D pdbqts with RDKit
    with tqdm_joblib(tqdm(desc="Generating...", total=len(smi_dict))) as progress_bar:
        Parallel(n_jobs=params.threads)(delayed(make_pdbs_from_smiles)(smi) for smi in smi_dict.items())

    # get a list of all successfully converted pdbs
    pdbs = os.listdir(os.path.join('utils','temp','pdb_files',''))

    logging.info('Converting pdbs to pdbqts...')

    # parallelise conversion of pdbs to pdbqt files
    with tqdm_joblib(tqdm(desc="Converting...", total=len(pdbs))) as progress_bar:
        Parallel(n_jobs=params.threads)(delayed(autodock_convert)(pdb_file, os.path.join('utils','MGLTools-1.5.6','')) for pdb_file in pdbs)

    # get list of successfully converted pdbqts
    pdbqts = get_filepaths(os.path.join('utils','temp','pdbqt_files',''))

    # get os name
    if sys.platform.lower() == 'darwin':
        os_name = 'mac'
    elif 'linux' in sys.platform.lower():
        os_name = 'linux'

    logging.info("Docking pdbqt ligands...")

    # then for each pdbqt
    for pdbqt in tqdm(pdbqts):

        # then dock each one
        dock_file(
                    os.path.join('utils','gwovina-1.0','build',os_name,'release','gwovina'),
                    params.receptor,
                    pdbqt,
                    *coords,
                    dock_settings['gwovina_settings']['exhaustiveness'],
                    dock_settings['gwovina_settings']['num_wolves'],
                    dock_settings['gwovina_settings']['num_modes'],
                    dock_settings['gwovina_settings']['energy_range'],
                    outfile=os.path.join(f'{stem_path}','utils','temp','docked_pdbqt_files',f'{os.path.split(pdbqt)[1]}')
                    )

    # get name of the input smiles file
    if '.' in params.ligand:
        docked_ligands_folder = os.path.basename(params.ligand).split('.')[0]
    else:
        docked_ligands_folder = os.path.basename(params.ligand)

    # define folder with the input file name to store the docked ligands
    docked_ligands_path = os.path.join('docked_ligands',docked_ligands_folder,'')

    # add each docked ligand in temp file as a ligand to score to the list of ligands in params.ligand
    params.ligand = [os.path.join('utils','temp','docked_pdbqt_files', file) for file in os.listdir(os.path.join('utils','temp','docked_pdbqt_files'))]
    
    # build receptor as a repeating list into params dict
    receptors = [params.receptor for i in range(len(params.ligand))]
    params.receptor = receptors

    # make sure the docked ligands folder exists
    if not os.path.isdir('docked_ligands'):
        os.mkdir('docked_ligands')
    if not os.path.isdir(docked_ligands_path):
        os.makedirs(docked_ligands_path)
    
    # copy the docked ligands into a main folder so they are accessible
    for file in params.ligand:
        shutil.copy(file, docked_ligands_path)   
    
    # return the updated parameters with docked ligands to score and the
    # dictionary of ids to smiles strings
    return params, smi_dict

def count_input_poses(list_of_ligands):

    """
    Function: Count the total number of poses across
              all ligands that need scoring

    Inputs:   List of ligand pdbqt files

    Outputs:  Total poses to score as integer
    """

    # set up count
    total_poses = 0

    # then for each ligand
    for ligand in tqdm(list_of_ligands):
        
        # open it and count the number of models
        ligand_file = open(ligand).read()
        poses = ligand_file.count('MODEL')

        # add the poses to total pose count
        if poses == 0:
            total_poses += 1
        else:
            total_poses += poses

    # return integer of total poses
    return total_poses

def ligand_pose_generator(params, lower_index, upper_index):

    """
    Function: Generate a list of receptor and ligand arguments
              between specific indices to pass to a function to be scored.
              This is necessary because if we do this at the start, if
              many ligands are being scored we cannot store them all in 
              memory and the script crashes

    Inputs:   Command line arguments, the lower and upper indexes
              of ligand poses to score

    Outputs:  List of tuples to score in the form
              [(receptor filename, ligand filename, (ligand pose number, ligand pdbqt block))]
    """

    # set up list to populate
    requested_receptor_ligand_args = list()

    # track the pose index
    pose_index = 0

    # then for each ligand
    for ligand_index, ligand_filepath in enumerate(params.ligand):

         
         # don't waste any time if we're already over the upper index
        if pose_index > upper_index:
            break

        # load the poses from the current ligand being considered
        pdbqt_pose_blocks = list()
        lig_text = open(ligand_filepath, 'r').read()
        lig_poses = lig_text.split('MODEL')
        for pose in lig_poses:
            lines = pose.split('\n')
            clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
            if len(clean_lines) < 3:
                pass
            else:
                pose = '\n'.join(clean_lines)
                pdbqt_pose_blocks.append(pose)

                # stop if we only want one pose
                if params.pose_1:
                    break

        # make a tuple with pdbqt block and pose name
        poses = [(f'_pose_{pdbqt_pose_blocks.index(pose) + 1}', pose) for pose in pdbqt_pose_blocks]

        # for each pose
        for pose in poses:
            
            # if the pose is one we want
            if lower_index <= pose_index < upper_index:
            
                # add it to the receptor ligand arguments
                receptor_ligand_args = (params.receptor[ligand_index], ligand_filepath, pose)

                requested_receptor_ligand_args.append(receptor_ligand_args)
    
            # update the pose index
            pose_index += 1

    # return the requested poses between the specified indexes
    return requested_receptor_ligand_args

def calculate_batches_needed(total_poses):

    """
    Function: Calculate how many batches the ligands
              need to be scored in with the available
              memory
    
    Input:    Total number of poses as an integer

    Output:   Number of batches to split the ligand poses into
    """

    # estimate the ram usage from empirical data
    estimated_ram_usage = (360540*total_poses) + 644792975
    available_ram = psutil.virtual_memory().total
    safe_ram_available = available_ram*0.7

    # use this info to split the batches into least possible
    # batches without getting a memory error
    if estimated_ram_usage > safe_ram_available:
        batches_needed = math.ceil(estimated_ram_usage/safe_ram_available)
    else:
        batches_needed = 1

    # return batch number as integer
    return batches_needed


def list_to_chunk_indexes(list_length, number_of_chunks):

    """
    Function: Create nested list of upper and lower indexes
              which relate to batches of input list

    Inputs:   List length (integer) and number of batches

    Outputs:  Nested list of chained indexes e.g.
              [(0, 573),(573, 1092)]
    """

    # set up list to populate
    indices = list()

    # get the size of each batch
    chunksize = math.ceil(list_length/number_of_chunks)

    # then for each chunk taking chunksize steps
    for i in range(0, list_length, chunksize):

        # if its the last chunk then it ends at the end of the list
        if i+chunksize < list_length:
            indices.append((i, i+chunksize))

        # otherwise it spans one chunksize chunk
        else:
            indices.append((i, list_length))

    # return list of indices
    return indices