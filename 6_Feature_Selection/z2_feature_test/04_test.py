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
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# os.environ['NUMEXPR_MAX_THREADS'] = '1'

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




def run_binana(pose_info,receptor_file):

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
    main_binana_out = binana.Binana(pose_info,receptor_file).out
    print(main_binana_out)


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

    #1. add closest contacts to binana_features dict(2.5)
    for contact in keep_closest_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['closest'].get(binana_name)
    
    #2. add close contacts to binana_features dict(4)
    for contact in keep_close_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['close'].get(binana_name)
    
    #3. add ligand atoms to binana_features dict as binary tallies
    for atom in keep_ligand_atoms:
        binana_name = atom.split()[-1]
        if main_binana_out['ligand_atoms'].get(binana_name) is None:
            binana_features[atom] = 0
        else:
            binana_features[atom] = 1
    
    #4. add electrostatics to binana_features dict
    for elsum in keep_elsums:
        binana_name = elsum.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[elsum] = main_binana_out['elsums'].get(binana_name)

    #6. add active site flexibility features to binana_features
    binana_features["BPF ALPHA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_ALPHA")
    binana_features["BPF ALPHA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_ALPHA")
    binana_features["BPF BETA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_BETA")
    binana_features["BPF BETA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_BETA")
    binana_features["BPF OTHER SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_OTHER")
    binana_features["BPF OTHER BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_OTHER")

    #8. add hydrophobic features to binana_features
    binana_features["HC ALPHA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_ALPHA")
    binana_features["HC ALPHA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_ALPHA")
    binana_features["HC BETA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_BETA")
    binana_features["HC BETA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_BETA")
    binana_features["HC OTHER SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_OTHER")
    binana_features["HC OTHER BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_OTHER")

    #7. add hydrogen bond features to binana_features
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

    #12. add salt bridge features to binana_features
    binana_features["SB ALPHA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_ALPHA")
    binana_features["SB BETA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_BETA")
    binana_features["SB OTHER"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_OTHER")

    #9.10 add aromatic stacking features to binana_features
    binana_features["piStack ALPHA"] = main_binana_out['stacking'].get("STACKING ALPHA")
    binana_features["piStack BETA"] = main_binana_out['stacking'].get("STACKING BETA")
    binana_features["piStack OTHER"] = main_binana_out['stacking'].get("STACKING OTHER")
    binana_features["tStack ALPHA"] = main_binana_out['t_stacking'].get("T-SHAPED_ALPHA")
    binana_features["tStack BETA"] = main_binana_out['t_stacking'].get("T-SHAPED_BETA")
    binana_features["tStack OTHER"] = main_binana_out['t_stacking'].get("T-SHAPED_OTHER")

    #11. add cation pi features to binana_features
    binana_features["catPi BETA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binana_features["catPi OTHER LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")

    #5. add rotatable bond count to binana features
    binana_features["nRot"] = main_binana_out['nrot']

    # return dictionary
    return binana_features



# 从pdbqt文件中读取配体姿态信息
ligand_file= "/home/s2331261/Master_Project/6_Feature_Selection/C1_get_decoy_poses/01_decoy_split_poses/1a0q_decoy_231_pose_3.pdbqt"
receptor_folder = "/home/s2331261/Master_Project/z1_3p_dataset/dock_source/all_receptors_ligands"
filename = os.path.basename(ligand_file)
pdb_code = filename.split('_')[0]
receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")


lig_text = open(ligand_file, 'r').read()
lines = lig_text.split('\n')
clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
if len(clean_lines) < 3:
    pass
else:
    pose_info = '\n'.join(clean_lines)


binana_dict = run_binana(pose_info,receptor_file)
binana_df = pd.DataFrame([binana_dict])
print("#######################################################################")
print(binana_df)
binana_df.to_csv('test_04_1.csv', index=False)
