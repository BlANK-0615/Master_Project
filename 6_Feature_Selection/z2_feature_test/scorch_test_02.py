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

# @contextlib.contextmanager
# def tqdm_joblib(tqdm_object):

#     """
#     Context manager to allow joblib
#     progress monitoring
#     """

#     def tqdm_print_progress(self):
#         if self.n_completed_tasks > tqdm_object.n:
#             n_completed = self.n_completed_tasks - tqdm_object.n
#             tqdm_object.update(n=n_completed)

#     original_print_progress = parallel.Parallel.print_progress
#     parallel.Parallel.print_progress = tqdm_print_progress

#     try:
#         yield tqdm_object
#     finally:
#         parallel.Parallel.print_progress = original_print_progress
#         tqdm_object.close()


# 从pdbqt文件中读取配体姿态信息
lig_text = open(file, 'r').read()
lines = lig_text.split('\n')
clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
if len(clean_lines) < 3:
    pass
else:
    pose_info = '\n'.join(clean_lines)



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
            

                ################################################################################33
                #REALLY IMPORTANT！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！1            
                # add it to the receptor ligand arguments
                receptor_ligand_args = (params.receptor[ligand_index], ligand_filepath, pose)

                requested_receptor_ligand_args.append(receptor_ligand_args)
    
                ################################################################################33
            # update the pose index
            pose_index += 1

    # return the requested poses between the specified indexes
    return requested_receptor_ligand_args


def prepare_features(receptor_ligand_args):

    filterwarnings('ignore')

    """
    Function: Wrapper to prepare all requested protein-ligand            
              complexes/poses for scoring             
                                            
    Inputs:   A tuple of a receptor ligand pair to score                            
                                            
    Output:   A single row dataframe of protein-ligand complex
              features                      
    """

    # grab the ligand pose number and pdbqt string block
    ligand_pose_number = receptor_ligand_args[2][0]
    ligand_pdbqt_block = receptor_ligand_args[2][1]
     
    # grab the ligand and receptor filepaths and filenames
    receptor_filepath = receptor_ligand_args[0]
    ligand_filepath = receptor_ligand_args[1]
    ligand_basename = os.path.basename(ligand_filepath)
    ligand_basename = ligand_basename.replace('.pdbqt', ligand_pose_number)
    receptor_basename = os.path.basename(receptor_filepath)

    # extract the interaction features
    features = extract(ligand_pdbqt_block, receptor_filepath)

    # prune the headers down to those needed for model scoring
    multi_pose_features = prune_df_headers(features)

    # fill None values with 0 for binana features
    multi_pose_features.fillna(0, inplace=True)

    # add receptor and ligand info to features
    multi_pose_features['Receptor'] = receptor_basename
    multi_pose_features['Ligand'] = ligand_basename
    
    # return the dataframe
    return multi_pose_features


# def score_ligand_batch(params, ligand_batch, model_binaries):

#     """
#     Function:  Make a dataframe of scores and stats for a batch of
#                ligands

#     Inputs:    Parsed command line parameters, batch of receptor ligand tuples,
#                and the binary files for the models to use
    
#     Outputs:   Dataframe of scores and stats for ligands in the batch
#     """
    ###############################################################
    # multiprocess the extracting of features from the protein ligand pairs
    # ligand_batch是个元组，每个元组里面有三个元素，分别是receptor文件名，ligand文件名，ligand的pose编号和pdbqt块
    with tqdm_joblib(tqdm(desc="Preparing features", total=len(ligand_batch))) as progress_bar:
        multi_pose_features = Parallel(n_jobs=params.threads)(delayed(prepare_features)(ligand) for ligand in ligand_batch)
    #####################################################################



# def scoring(params):

#     """
#     Function: Score protein-ligand complex(es)                             
                                            
#     Inputs: User command line parameters dictionary                              
                                            
#     Output: Dataframe of scoring function predictions                             
#     """

#     # print help/intro if requested
#     print_intro(params)

#     # prepare and dock smiles if smiles ligands supplied
#     if params.dock:

        ##############################################################################################################
        # load in ligands to score for this batch
        ligand_batch = ligand_pose_generator(params, index_range[0], index_range[1])
        #ligand_batch是个元组，每个元组里面有三个元素，分别是receptor文件名，ligand文件名，ligand的pose编号和pdbqt块
        ##############################################################################################################
        
#         # score the batch and get the results
#         merged_results = score_ligand_batch(params, ligand_batch, model_binaries)

#         # add the results to the ligand_scores list
#         ligand_scores.append(merged_results)

#     # make final results from ligand scores
#     final_ligand_scores = create_final_results(params, ligand_scores)

#     # if we have docked the ligands then add column for their smiles strings
#     if params.dock:
#         final_ligand_scores['Ligand_SMILE'] = final_ligand_scores['Ligand_ID'].map(smi_dict)


#     # return the final dataframe of ligand scores
#     return final_ligand_scores

#######################################################################
# Main script

if __name__ == "__main__":

    # parse user arguments
    params = parse_args(sys.argv)

    # score complexes
    scoring_function_results = scoring(params)

    # output results
    if not params.out:

        # send to stdout if no outfile given
        sys.stdout.write(scoring_function_results.to_csv(index=False))

    else:

        # otherwise save to user specified csv
        scoring_function_results.to_csv(params.out, index=False)

# scorch.py end
# 明天要看下binana的输出是什么，如果可以直接做，就按这个脚本直接做，直接合并成一个dataframe