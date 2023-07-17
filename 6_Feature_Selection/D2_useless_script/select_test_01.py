"""
SCORCH script version 1.0 - Run python scoring.py -h for help    
                                                                    
Script Authors:                                                     
@sammoneykyrle                                                      
@milesmcgibbon                                                      
                                                                    
School of Biological Sciences                                       
The University of Edinburgh                                         
"""

#######################################################################
# Imports

# import os and set tensorflow verbosity
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


def run_binana(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get BINANA descriptors for    
    protein-ligand complex                  
                                        
    Inputs:  ligand as a pdbqt string block,         
    receptor pdbqt filepath                 
                                            
    Output: BINANA protein-ligand complex   
    descriptor features as a dictionary     
    """

    # empty dictionary to populate with features
    binana_features = dict()

    # get dictionary of binana features
    main_binana_out = binana.Binana(ligand_pdbqt_block, receptor_filepath).out


    # define the features we want
    keep_closest_contacts = ["2.5 (HD, OA)", 
                             "2.5 (HD, HD)", 
                             "2.5 (HD, N)", 
                             "2.5 (C, HD)", 
                             "2.5 (OA, ZN)", 
                             "2.5 (HD, ZN)", 
                             "2.5 (A, HD)"]
    
    keep_close_contacts = ["4.0 (C, C)", 
                           "4.0 (HD, OA)", 
                           "4.0 (C, HD)", 
                           "4.0 (C, N)", 
                           "4.0 (A, C)",
                           "4.0 (A, OA)", 
                           "4.0 (N, OA)", 
                           "4.0 (A, N)", 
                           "4.0 (HD, N)", 
                           "4.0 (HD, HD)", 
                           "4.0 (A, HD)", 
                           "4.0 (OA, OA)", 
                           "4.0 (C, OA)", 
                           "4.0 (N, N)",
                           "4.0 (C, SA)", 
                           "4.0 (HD, SA)", 
                           "4.0 (OA, SA)", 
                           "4.0 (N, SA)", 
                           "4.0 (A, A)", 
                           "4.0 (HD, S)", 
                           "4.0 (S, ZN)", 
                           "4.0 (N, ZN)", 
                           "4.0 (HD, ZN)", 
                           "4.0 (A, SA)", 
                           "4.0 (OA, ZN)", 
                           "4.0 (C, ZN)", 
                           "4.0 (C, NA)", 
                           "4.0 (NA, OA)", 
                           "4.0 (HD, NA)", 
                           "4.0 (N, NA)", 
                           "4.0 (A, NA)", 
                           "4.0 (BR, C)", 
                           "4.0 (HD, P)", 
                           "4.0 (F, N)", 
                           "4.0 (F, HD)", 
                           "4.0 (C, CL)", 
                           "4.0 (CL, HD)"]

    keep_ligand_atoms = ["LA N",
                         "LA HD"]
    
    keep_elsums = [ "ElSum (C, C)",
                    "ElSum (HD, OA)",
                    "ElSum (C, HD)",
                    "ElSum (C, N)",
                    "ElSum (A, C)",
                    "ElSum (A, OA)",
                    "ElSum (N, OA)",
                    "ElSum (A, N)",
                    "ElSum (HD, HD)",
                    "ElSum (A, HD)",
                    "ElSum (OA, OA)",
                    "ElSum (C, OA)",
                    "ElSum (N, N)",
                    "ElSum (C, SA)",
                    "ElSum (HD, SA)",
                    "ElSum (OA, SA)",
                    "ElSum (N, SA)",
                    "ElSum (A, A)",
                    "ElSum (N, S)",
                    "ElSum (HD, S)",
                    "ElSum (OA, S)",
                    "ElSum (A, SA)",
                    "ElSum (C, NA)",
                    "ElSum (NA, OA)",
                    "ElSum (HD, NA)",
                    "ElSum (N, NA)",
                    "ElSum (A, NA)",
                    "ElSum (BR, C)",
                    "ElSum (HD, P)",
                    "ElSum (OA, P)",
                    "ElSum (N, P)",
                    "ElSum (C, F)",
                    "ElSum (F, N)",
                    "ElSum (A, F)",
                    "ElSum (CL, OA)",
                    "ElSum (C, CL)",
                    "ElSum (CL, N)",
                    "ElSum (A, CL)"]

    # add closest contacts to binana_features dict
    for contact in keep_closest_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['closest'].get(binana_name)
    
    # add close contacts to binana_features dict
    for contact in keep_close_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['close'].get(binana_name)
    
    # add ligand atoms to binana_features dict as binary tallies
    for atom in keep_ligand_atoms:
        binana_name = atom.split()[-1]
        if main_binana_out['ligand_atoms'].get(binana_name) is None:
            binana_features[atom] = 0
        else:
            binana_features[atom] = 1
    
    # add electrostatics to binana_features dict
    for elsum in keep_elsums:
        binana_name = elsum.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[elsum] = main_binana_out['elsums'].get(binana_name)

    # add active site flexibility features to binana_features
    binana_features["BPF ALPHA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_ALPHA")
    binana_features["BPF ALPHA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_ALPHA")
    binana_features["BPF BETA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_BETA")
    binana_features["BPF BETA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_BETA")
    binana_features["BPF OTHER SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_OTHER")
    binana_features["BPF OTHER BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_OTHER")

    # add hydrophobic features to binana_features
    binana_features["HC ALPHA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_ALPHA")
    binana_features["HC ALPHA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_ALPHA")
    binana_features["HC BETA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_BETA")
    binana_features["HC BETA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_BETA")
    binana_features["HC OTHER SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_OTHER")
    binana_features["HC OTHER BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_OTHER")

    # add hydrogen bond features to binana_features
    binana_features["HB ALPHA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_ALPHA")
    binana_features["HB BETA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_OTHER")
    binana_features["HB ALPHA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_ALPHA")
    binana_features["HB ALPHA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_ALPHA")
    binana_features["HB BETA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_OTHER")

    # add salt bridge features to binana_features
    binana_features["SB ALPHA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_ALPHA")
    binana_features["SB BETA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_BETA")
    binana_features["SB OTHER"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_OTHER")

    # add aromatic stacking features to binana_features
    binana_features["piStack ALPHA"] = main_binana_out['stacking'].get("STACKING ALPHA")
    binana_features["piStack BETA"] = main_binana_out['stacking'].get("STACKING BETA")
    binana_features["piStack OTHER"] = main_binana_out['stacking'].get("STACKING OTHER")
    binana_features["tStack ALPHA"] = main_binana_out['t_stacking'].get("T-SHAPED_ALPHA")
    binana_features["tStack BETA"] = main_binana_out['t_stacking'].get("T-SHAPED_BETA")
    binana_features["tStack OTHER"] = main_binana_out['t_stacking'].get("T-SHAPED_OTHER")

    # add cation pi features to binana_features
    binana_features["catPi BETA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binana_features["catPi OTHER LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")

    # add rotatable bond count to binana features
    binana_features["nRot"] = main_binana_out['nrot']

    # return dictionary
    return binana_features

def kier_flexibility(ligand_pdbqt_block):

    """
    Function: Calculate Kier flexibility    
    for ligand                              
                                            
    Inputs: ligand as a pdbqt string block  
                                            
    Output: Kier flexibility                
    """

    # parse pdbqt block
    mol = kier.SmilePrep(ligand_pdbqt_block)

    # calculate flexibility
    return kier.CalculateFlexibility(mol)

def calculate_ecifs(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get ECIFs for protein-ligand  
    complex                                 
                                            
    Inputs: ligand as a pdbqt string block, 
    receptor pdbqt filepath                 
                                            
    Output: ECIF protein-ligand complex     
    descriptor features as a pandas DataFrame      
    """

    # get ECIFs with default cutoff using imported functions
    # from utils/ecifs.py
    ECIF_data = GetECIF(receptor_filepath, ligand_pdbqt_block, distance_cutoff=6.0)

    # replace the semicolons to make valid dataframe headers
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]

    # zip into a dictionary and convert to dataframe
    ECIF_data = dict(zip(ECIFHeaders,ECIF_data))
    ECIF_df = pd.DataFrame(ECIF_data,index=[0])

    # return the dataframe
    return ECIF_df

def extract(ligand_pdbqt_block, receptor_filepath):

    """
    Function: Get all descriptor features   
    for protein-ligand complex              
                                            
    Inputs: ligand as a pdbqt string block, 
    receptor pdbqt filepath     
                                            
    Output: All protein-ligand complex      
    descriptor features as a DataFrame      
    """
    
    # get the kier flexibility
    k = kier_flexibility(ligand_pdbqt_block)

    # get the binana descriptors and build into dataframe
    binana_dict = run_binana(ligand_pdbqt_block,receptor_filepath)
    binana_df = pd.DataFrame([binana_dict])

    # get the ECIFs
    ECIF = calculate_ecifs(ligand_pdbqt_block, receptor_filepath)

    # concatenate all feature columns to one row dataframe
    df = pd.concat([ECIF,binana_df],axis=1)

    # add the kier flexibility as a column
    df['Kier Flexibility'] = k

    # return the features
    return df

def prune_df_headers(df):

    """
    Function: Condense features for model input                             
                                            
    Inputs: Full Dataframe of               
    protein-ligand complex descriptors     
                                            
    Output: DataFrame of features for model input                                   
    """

    # load features we want in the correct order
    reference_headers = json.load(open(os.path.join('utils','params','features.json')))
    headers_58 = reference_headers.get('492_models_58')

    # subset the dataframe of features
    df = df[headers_58]

    # return the dataframe
    return df

def multiple_pose_check(ligand_filepath):

    """
    Function: Transform ligand.pdbqt file of      
    poses/models into pdbqt string blocks   
                                            
    Inputs: ligand.pdbqt filepath           
                                            
    Output: List of model/pose pdbqt string 
    blocks as a tuple with the pose number
    e.g. [('_pose_1','REMARK....)]                                 
    """

    # make empty list for populating
    pdbqt_pose_blocks = list()

    # open the input ligand file
    lig_text = open(ligand_filepath, 'r').read()

    # split by poses
    lig_poses = lig_text.split('MODEL')

    # for each pose clean up any whitespace or empty lines
    for pose in lig_poses:
        lines = pose.split('\n')
        clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]

        # if there are less than three lines then
        # its an artefact of the splitting and we ignore 
        if len(clean_lines) < 3:
            pass

        # otherwise join up the correct poses
        else:
            pose = '\n'.join(clean_lines)

            # add the pose to the list of poses
            pdbqt_pose_blocks.append(pose)

    # map the poses into a list of tuples
    pdbqt_pose_blocks = list(map(lambda x: (f'_pose_{pdbqt_pose_blocks.index(x) + 1}', x), pdbqt_pose_blocks))

    # return the list of tuples
    return pdbqt_pose_blocks

def run_networks(df, models_to_load, model_name):

    """
    Function: Get single mean prediction from consensus of 
    neural network models for a protein-ligand complex              
                                            
    Inputs: Dataframe of features for model,
            list of model filenames to load,
            model name string for prediction columns

    Output: Dataframe of model predictions with 
            final average mean column              
    """

    # make empty dataframe to populate wit predictions
    predictions = pd.DataFrame()

    # for each model
    for i in tqdm(range(len(models_to_load))):

        # load the model and get its prediction on the data
        model = load_model(models_to_load[i])
        y_pred = model.predict(df)

        # add the prediction as a column to the predictions df
        predictions[f'{model_name}_{i + 1}'] = y_pred.flatten()

    # take an average of all the predictions
    predictions[f'{model_name}_models_average'] = predictions.mean(axis=1)

    # return the df of predictions
    return predictions


def binary_concat(dfs, headers):

    """
    Function: Concatenate list of           
    dataframes into a single dataframe by   
    sequentially writing to a single binary 
    file (removes pd.concat bottleneck)     
                                            
    Inputs: List of dataframes, dataframe   
    headers as a list                       
                                            
    Output: Single combined dataframe       
    """

    # set total rows count
    total_rows = 0

    # make temporary directory for binary file if it
    # doesn't exist yet
    if not os.path.isdir(os.path.join('utils','temp')):
        os.makedirs(os.path.join('utils','temp'))

    # create a temporary binary file
    with open(os.path.join('utils','temp','features.bin'),'wb+') as binary_store:

        # then for each dataframe write it to binary
        for df in dfs:

            # make sure nRot is numeric in type
            df['nRot'] = pd.to_numeric(df['nRot'])

            # get the rows and columns of the dataframe
            rows, fixed_total_columns = df.shape

            # add the rows to total rows counter
            total_rows += rows

            # write the df to the binary store file
            binary_store.write(df.values.tobytes())

            # store the dtypes
            typ = df.values.dtype

    # open the temporary binary file
    with open(os.path.join('utils','temp','features.bin'),'rb') as binary_store:

        # read the binary file
        buffer = binary_store.read()

        # shape the data into a dataframe using the info from
        # when we wrote all the files to the binary store
        data = np.frombuffer(buffer, dtype=typ).reshape(total_rows, fixed_total_columns)
        master_df = pd.DataFrame(data = data, columns = headers)

    # make sure the binary store is deleted in case we call the function again
    os.remove(os.path.join('utils','temp','features.bin'))

    # return the concatenated dataframe
    return master_df

def parse_module_args(args_dict):

    """
    Function: Parse user arguments when     
    script is imported as a module          
                                            
    Inputs: User arguments as a dictionary  
                                            
    Output: Populated params dictionary     
    """

    # empty list to use as a spoof sys.argv result
    command_input = list()

    # check if any boolean flag arguments have been passed
    boolean_args = ['verbose','return_pose_scores']
    for key, value in args_dict.items():
        if key in boolean_args:
            if value:
                command_input.append(f'-{key}')
        
        # otherwise add them as normal args with their values
        else:
            command_input.append(f'-{key}')
            command_input.append(str(value))

    # then parse the args as if they were command line
    parsed_args = parse_args(command_input)

    # return the arguments
    return parsed_args

def parse_args(args):

    """
    Function: Parse user defined command    
    line arguments                          
                                            
    Inputs: Command line arguments          
                                            
    Output: Populated params dictionary     
    """

    parser = argparse.ArgumentParser(description="SCORCH 1.0\nMiles McGibbon, Samuel Money-Kyrle, Vincent Blay & Douglas R. Houston")

    requiredNamed = parser.add_argument_group('required named arguments')
     
    # add required arguments
    requiredNamed.add_argument('-l','--ligand', help="""Ligands to score against the supplied receptor. 
                                                 Can be a .smi or .txt filepath, a .pdbqt filepath, or path to a folder of pdbqt files. 
                                                 If .smi file is supplied, --range and --center args or --ref_lig args are 
                                                 also required.""", required=True)
    requiredNamed.add_argument('-r','--receptor', help="Receptor to score ligands against. Must be a filepath to a .pdbqt file", required=True)

    # add optional arguments
    parser.add_argument('-rl','--ref_lig', help="Filepath to example ligand in receptor binding site (mol, mol2, sdf, pdb or pdbqt)")
    parser.add_argument('-t','--threads', default=1, help="Number of CPU threads to parallelise SCORCH over", type=int)
    parser.add_argument('-c','--center', help="'[x, y, z]' coordinates of the center of the binding site for docking")
    parser.add_argument('-ra','--range', help="'[x, y, z]' axis lengths to define a box around --center coordinates for docking")
    parser.add_argument('-o','--out', help="Filepath for output csv (If not supplied, scores are written to stdout)")
    parser.add_argument('-p','--return_pose_scores', action='store_true', help="If supplied, scoring values for individual poses in each ligand file are returned")
    parser.add_argument('-v','--verbose', action='store_true', help="If supplied, progress bars and indicators are displayed while scoring")
    parser.add_argument('-s','--pose_1', action='store_true', help="Consider only the first pose in each pdbqt file to score - NOT RECOMMENDED")
    parser.add_argument('-d','--dock', action='store_true', help="""If supplied, input ligands are assumed to be text SMILES and will be docked 
                                                                    using GWOVina before scoring. This will be autodetected if a .smi or .txt file is supplied""")
    params = parser.parse_args()

    if params.verbose:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    else:
        tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)
        logging.basicConfig(level=logging.WARNING, format='%(message)s')

    if os.path.isdir(params.ligand):
        params.ligand = [os.path.join(params.ligand, file) for file in os.listdir(params.ligand)]
        receptors = [params.receptor for i in range(len(params.ligand))]
        params.receptor = receptors

    elif '.smi' in params.ligand or '.txt' in params.ligand:
        params.dock = True
    else:
        params.ligand = [params.ligand]
        params.receptor = [params.receptor]

    return params