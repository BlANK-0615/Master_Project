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
    keep_closest_contacts = ['2.5 (HD, OA)', '2.5 (HD, HD)', '2.5 (C, HD)', '2.5 (HD, NA)', '2.5 (FE, HD)', '2.5 (A, HD)', '2.5 (HD, N)', '2.5 (OA, ZN)', '2.5 (N, ZN)', '2.5 (HD, ZN)', '2.5 (OA, OA)', '2.5 (CO, OA)', '2.5 (CO, HD)', '2.5 (C, C)', '2.5 (N, OA)', '2.5 (F, HD)', '2.5 (C, OA)', '2.5 (C, CL)', '2.5 (MG, OA)', '2.5 (NA, OA)', '2.5 (CA, OA)', '2.5 (CL, HD)', '2.5 (CL, MG)', '2.5 (A, OA)', '2.5 (C, MG)', '2.5 (C, N)', '2.5 (HD, SA)', '2.5 (MN, OA)', '2.5 (FE, NA)', '2.5 (FE, OA)', '2.5 (HD, P)', '2.5 (CA, HD)', '2.5 (HD, MN)', '2.5 (A, P)', '2.5 (F, OA)', '2.5 (C, ZN)', '2.5 (HD, MG)', '2.5 (N, N)', '2.5 (C, NA)', '2.5 (CO, NA)', '2.5 (A, C)', '2.5 (NA, ZN)', '2.5 (A, CL)', '2.5 (MG, N)', '2.5 (CA, NA)', '2.5 (MN, N)', '2.5 (CD, OA)', '2.5 (A, N)', '2.5 (CL, SA)', '2.5 (OA, SA)', '2.5 (HD, S)', '2.5 (BR, HD)', '2.5 (SA, SA)', '2.5 (SA, ZN)', '2.5 (N, SA)', '2.5 (NI, OA)', '2.5 (C, CA)', '2.5 (A, SA)', '2.5 (C, SA)', '2.5 (N, NA)', '2.5 (CL, N)', '2.5 (C, P)', '2.5 (A, A)', '2.5 (MN, NA)', '2.5 (CU, OA)', '2.5 (HD, NI)', '2.5 (C, CO)', '2.5 (A, MG)', '2.5 (CU, HD)', '2.5 (CO, N)', '2.5 (A, NA)', '2.5 (C, FE)', '2.5 (CL, OA)', '2.5 (N, NI)', '2.5 (BR, OA)', '2.5 (BR, C)', '2.5 (A, ZN)', '2.5 (C, F)', '2.5 (MG, P)', '2.5 (OA, P)', '2.5 (A, F)', '2.5 (MG, NA)', '2.5 (C, NI)', '2.5 (HD, K)', '2.5 (BR, FE)', '2.5 (A, BR)', '2.5 (C, MN)', '2.5 (F, N)', '2.5 (CU, NA)', '2.5 (CD, HD)', '2.5 (CU, N)', '2.5 (A, MN)', '2.5 (A, FE)', '2.5 (NA, NA)', '2.5 (F, MG)', '2.5 (A, CU)', '2.5 (CU, F)', '2.5 (FE, N)', '2.5 (C, CU)', '2.5 (NA, NI)', '2.5 (C, CD)', '2.5 (F, MN)', '2.5 (K, OA)', '2.5 (C, K)', '2.5 (S, ZN)', '2.5 (K, NA)', '2.5 (F, ZN)', '2.5 (CD, NA)', '2.5 (N, P)', '2.5 (F, SA)', '2.5 (F, FE)', '2.5 (OA, S)', '2.5 (BR, SA)', '2.5 (MG, SA)', '2.5 (BR, ZN)', '2.5 (A, NI)', '2.5 (CL, ZN)', '2.5 (CA, N)', '2.5 (MN, P)', '2.5 (HD, I)', '2.5 (P, ZN)', '2.5 (NA, SA)']
    #len(keep_closest_contacts)
    #122
    keep_close_contacts = ['4.0 (C, C)', '4.0 (C, OA)', '4.0 (C, HD)', '4.0 (C, N)', '4.0 (N, OA)', '4.0 (HD, OA)', '4.0 (N, N)', '4.0 (HD, N)', '4.0 (HD, HD)', '4.0 (A, OA)', '4.0 (A, C)', '4.0 (A, A)', '4.0 (A, N)', '4.0 (A, HD)', '4.0 (C, SA)', '4.0 (HD, SA)', '4.0 (OA, OA)', '4.0 (C, ZN)', '4.0 (OA, ZN)', '4.0 (C, NA)', '4.0 (A, NA)', '4.0 (NA, OA)', '4.0 (N, NA)', '4.0 (HD, NA)', '4.0 (A, SA)', '4.0 (OA, SA)', '4.0 (BR, C)', '4.0 (BR, OA)', '4.0 (F, OA)', '4.0 (F, HD)', '4.0 (A, F)', '4.0 (F, N)', '4.0 (C, F)', '4.0 (A, FE)', '4.0 (C, FE)', '4.0 (FE, N)', '4.0 (FE, NA)', '4.0 (FE, OA)', '4.0 (C, P)', '4.0 (CL, N)', '4.0 (CL, HD)', '4.0 (C, CL)', '4.0 (CL, OA)', '4.0 (N, ZN)', '4.0 (HD, ZN)', '4.0 (N, SA)', '4.0 (SA, ZN)', '4.0 (OA, P)', '4.0 (N, P)', '4.0 (HD, P)', '4.0 (MG, OA)', '4.0 (A, ZN)', '4.0 (HD, S)', '4.0 (A, S)', '4.0 (N, S)', '4.0 (OA, S)', '4.0 (S, ZN)', '4.0 (A, P)', '4.0 (A, CL)', '4.0 (P, ZN)', '4.0 (C, S)', '4.0 (CO, P)', '4.0 (CO, HD)', '4.0 (F, ZN)', '4.0 (F, SA)', '4.0 (S, SA)', '4.0 (FE, HD)', '4.0 (NA, ZN)', '4.0 (CO, N)', '4.0 (A, CO)', '4.0 (C, CO)', '4.0 (CO, OA)', '4.0 (HD, MG)', '4.0 (C, MG)', '4.0 (CL, SA)', '4.0 (NA, SA)', '4.0 (SA, SA)', '4.0 (A, CA)', '4.0 (C, CA)', '4.0 (CA, N)', '4.0 (CA, OA)', '4.0 (CA, P)', '4.0 (CA, HD)', '4.0 (A, MG)', '4.0 (C, MN)', '4.0 (MN, N)', '4.0 (HD, MN)', '4.0 (MN, NA)', '4.0 (A, MN)', '4.0 (F, MN)', '4.0 (BR, N)', '4.0 (BR, HD)', '4.0 (A, BR)', '4.0 (MN, P)', '4.0 (MN, OA)', '4.0 (C, I)', '4.0 (I, SA)', '4.0 (A, I)', '4.0 (MG, NA)', '4.0 (MG, N)', '4.0 (MG, P)', '4.0 (HD, NI)', '4.0 (A, NI)', '4.0 (NI, OA)', '4.0 (I, OA)', '4.0 (C, NI)', '4.0 (MG, SA)', '4.0 (N, NI)', '4.0 (C, CD)', '4.0 (CD, N)', '4.0 (F, MG)', '4.0 (F, S)', '4.0 (C, CU)', '4.0 (P, SA)', '4.0 (BR, SA)', '4.0 (CL, ZN)', '4.0 (MN, SA)', '4.0 (BR, MG)', '4.0 (CA, F)', '4.0 (FE, P)', '4.0 (CL, CU)', '4.0 (CA, NA)', '4.0 (NA, NI)', '4.0 (CU, N)', '4.0 (CU, HD)', '4.0 (CU, OA)', '4.0 (MG, S)', '4.0 (A, CU)', '4.0 (NA, P)', '4.0 (CO, NA)', '4.0 (CU, NA)', '4.0 (BR, ZN)', '4.0 (HD, I)', '4.0 (I, N)', '4.0 (F, NA)', '4.0 (NA, NA)', '4.0 (FE, SA)', '4.0 (C, K)', '4.0 (K, OA)', '4.0 (K, N)', '4.0 (HD, K)', '4.0 (K, P)', '4.0 (FE, S)', '4.0 (CD, OA)', '4.0 (CD, HD)', '4.0 (CA, SA)', '4.0 (A, HG)', '4.0 (HG, SA)', '4.0 (CD, P)', '4.0 (CO, S)', '4.0 (CA, S)', '4.0 (CL, MG)', '4.0 (CU, F)', '4.0 (NI, SA)', '4.0 (MN, S)', '4.0 (BR, NA)', '4.0 (F, FE)', '4.0 (F, K)', '4.0 (K, SA)', '4.0 (CL, FE)', '4.0 (A, CD)', '4.0 (CL, MN)', '4.0 (NA, S)', '4.0 (CA, CL)', '4.0 (BR, FE)', '4.0 (A, K)', '4.0 (K, NA)', '4.0 (CD, NA)', '4.0 (CL, NA)', '4.0 (CO, SA)', '4.0 (HG, OA)', '4.0 (CL, NI)', '4.0 (CD, SA)', '4.0 (CD, S)', '4.0 (BR, MN)']
    #len(keep_close_contacts)
    #175
    keep_ligand_atoms = ['LA C', 'LA N', 'LA HD', 'LA OA', 'LA A', 'LA SA', 'LA NA', 'LA BR', 'LA F', 'LA P', 'LA CL', 'LA S', 'LA I']
    #len(keep_ligand_atoms)
    #13
    keep_elsums = ['ElSum (C, C)', 'ElSum (C, OA)', 'ElSum (C, HD)', 'ElSum (C, N)', 'ElSum (N, OA)', 'ElSum (HD, OA)', 'ElSum (N, N)', 'ElSum (HD, N)', 'ElSum (HD, HD)', 'ElSum (A, OA)', 'ElSum (A, C)', 'ElSum (A, A)', 'ElSum (A, N)', 'ElSum (A, HD)', 'ElSum (C, SA)', 'ElSum (HD, SA)', 'ElSum (OA, OA)', 'ElSum (C, ZN)', 'ElSum (OA, ZN)', 'ElSum (C, NA)', 'ElSum (A, NA)', 'ElSum (NA, OA)', 'ElSum (N, NA)', 'ElSum (HD, NA)', 'ElSum (A, SA)', 'ElSum (OA, SA)', 'ElSum (BR, C)', 'ElSum (BR, OA)', 'ElSum (F, OA)', 'ElSum (F, HD)', 'ElSum (A, F)', 'ElSum (F, N)', 'ElSum (C, F)', 'ElSum (A, FE)', 'ElSum (C, FE)', 'ElSum (FE, N)', 'ElSum (FE, HD)', 'ElSum (FE, NA)', 'ElSum (FE, OA)', 'ElSum (C, P)', 'ElSum (CL, N)', 'ElSum (CL, HD)', 'ElSum (C, CL)', 'ElSum (CL, OA)', 'ElSum (N, ZN)', 'ElSum (HD, ZN)', 'ElSum (N, SA)', 'ElSum (SA, ZN)', 'ElSum (OA, P)', 'ElSum (N, P)', 'ElSum (HD, P)', 'ElSum (MG, OA)', 'ElSum (A, ZN)', 'ElSum (HD, S)', 'ElSum (A, S)', 'ElSum (N, S)', 'ElSum (OA, S)', 'ElSum (S, ZN)', 'ElSum (A, P)', 'ElSum (A, CL)', 'ElSum (P, ZN)', 'ElSum (C, S)', 'ElSum (CO, P)', 'ElSum (CO, OA)', 'ElSum (CO, HD)', 'ElSum (F, ZN)', 'ElSum (F, SA)', 'ElSum (S, SA)', 'ElSum (NA, ZN)', 'ElSum (CO, N)', 'ElSum (A, CO)', 'ElSum (C, CO)', 'ElSum (HD, MG)', 'ElSum (C, MG)', 'ElSum (CL, SA)', 'ElSum (NA, SA)', 'ElSum (SA, SA)', 'ElSum (A, CA)', 'ElSum (C, CA)', 'ElSum (CA, N)', 'ElSum (CA, OA)', 'ElSum (CA, P)', 'ElSum (CA, HD)', 'ElSum (A, MG)', 'ElSum (CL, MG)', 'ElSum (C, MN)', 'ElSum (MN, N)', 'ElSum (MN, OA)', 'ElSum (HD, MN)', 'ElSum (MN, NA)', 'ElSum (A, MN)', 'ElSum (F, MN)', 'ElSum (BR, N)', 'ElSum (BR, HD)', 'ElSum (A, BR)', 'ElSum (MN, P)', 'ElSum (C, I)', 'ElSum (I, SA)', 'ElSum (A, I)', 'ElSum (MG, NA)', 'ElSum (MG, N)', 'ElSum (MG, P)', 'ElSum (HD, NI)', 'ElSum (CO, NA)', 'ElSum (A, NI)', 'ElSum (NI, OA)', 'ElSum (I, OA)', 'ElSum (C, NI)', 'ElSum (MG, SA)', 'ElSum (CA, NA)', 'ElSum (N, NI)', 'ElSum (C, CD)', 'ElSum (CD, N)', 'ElSum (CD, OA)', 'ElSum (F, MG)', 'ElSum (F, S)', 'ElSum (C, CU)', 'ElSum (P, SA)', 'ElSum (BR, SA)', 'ElSum (CL, ZN)', 'ElSum (MN, SA)', 'ElSum (BR, MG)', 'ElSum (CA, F)', 'ElSum (FE, P)', 'ElSum (CL, CU)', 'ElSum (NA, NI)', 'ElSum (CU, N)', 'ElSum (CU, OA)', 'ElSum (CU, HD)', 'ElSum (MG, S)', 'ElSum (A, CU)', 'ElSum (NA, P)', 'ElSum (CU, NA)', 'ElSum (BR, ZN)', 'ElSum (HD, I)', 'ElSum (I, N)', 'ElSum (F, NA)', 'ElSum (NA, NA)', 'ElSum (FE, SA)', 'ElSum (C, K)', 'ElSum (K, OA)', 'ElSum (K, N)', 'ElSum (HD, K)', 'ElSum (K, P)', 'ElSum (FE, S)', 'ElSum (CD, HD)', 'ElSum (CA, SA)', 'ElSum (BR, FE)', 'ElSum (A, HG)', 'ElSum (HG, SA)', 'ElSum (CD, P)', 'ElSum (CO, S)', 'ElSum (CA, S)', 'ElSum (CU, F)', 'ElSum (NI, SA)', 'ElSum (MN, S)', 'ElSum (BR, NA)', 'ElSum (F, FE)', 'ElSum (F, K)', 'ElSum (K, SA)', 'ElSum (CL, FE)', 'ElSum (A, CD)', 'ElSum (CL, MN)', 'ElSum (NA, S)', 'ElSum (CA, CL)', 'ElSum (A, K)', 'ElSum (K, NA)', 'ElSum (CD, NA)', 'ElSum (CL, NA)', 'ElSum (CO, SA)', 'ElSum (HG, OA)', 'ElSum (CL, NI)', 'ElSum (CD, SA)', 'ElSum (CD, S)', 'ElSum (BR, MN)']
    #len(keep_elsums)
    #175
    #1.ok add closest contacts to binana_features dict(2.5)
    for contact in keep_closest_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['closest'].get(binana_name)
    
    #2.ok add close contacts to binana_features dict(4)
    for contact in keep_close_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['close'].get(binana_name)
    
    #3.ok add ligand atoms to binana_features dict as binary tallies
    for atom in keep_ligand_atoms:
        binana_name = atom.split()[-1]
        if main_binana_out['ligand_atoms'].get(binana_name) is None:
            binana_features[atom] = 0
        else:
            binana_features[atom] = 1
    
    #4.ok add electrostatics to binana_features dict
    for elsum in keep_elsums:
        binana_name = elsum.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[elsum] = main_binana_out['elsums'].get(binana_name)

    #6. ok add active site flexibility features to binana_features
    binana_features["BPF ALPHA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_ALPHA")
    binana_features["BPF ALPHA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_ALPHA")
    binana_features["BPF BETA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_BETA")
    binana_features["BPF BETA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_BETA")
    binana_features["BPF OTHER SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_OTHER")
    binana_features["BPF OTHER BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_OTHER")

    #8. ok add hydrophobic features to binana_features
    binana_features["HC ALPHA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_ALPHA")
    binana_features["HC ALPHA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_ALPHA")
    binana_features["HC BETA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_BETA")
    binana_features["HC BETA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_BETA")
    binana_features["HC OTHER SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_OTHER")
    binana_features["HC OTHER BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_OTHER")

    #7. ok add hydrogen bond features to binana_features
    binana_features["HB ALPHA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_ALPHA")
    binana_features["HB ALPHA BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_ALPHA")
    binana_features["HB BETA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_OTHER")

    #ok
    binana_features["HB ALPHA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_ALPHA")
    binana_features["HB ALPHA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_ALPHA")
    binana_features["HB BETA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_OTHER")

    #12. ok add salt bridge features to binana_features
    binana_features["SB ALPHA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_ALPHA")
    binana_features["SB BETA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_BETA")
    binana_features["SB OTHER"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_OTHER")

    #9.10 ok add aromatic stacking features to binana_features
    binana_features["piStack ALPHA"] = main_binana_out['stacking'].get("STACKING ALPHA")
    binana_features["piStack BETA"] = main_binana_out['stacking'].get("STACKING BETA")
    binana_features["piStack OTHER"] = main_binana_out['stacking'].get("STACKING OTHER")
    binana_features["tStack ALPHA"] = main_binana_out['t_stacking'].get("T-SHAPED_ALPHA")
    binana_features["tStack BETA"] = main_binana_out['t_stacking'].get("T-SHAPED_BETA")
    binana_features["tStack OTHER"] = main_binana_out['t_stacking'].get("T-SHAPED_OTHER")

    #11. add cation pi features to binana_features
    binana_features["catPi ALPHA RECEPTOR"] = main_binana_out['pi_cation'].get("PI-CATION_RECEPTOR-CHARGED_ALPHA")
    binana_features["catPi BETA RECEPTOR"] = main_binana_out['pi_cation'].get("PI-CATION_RECEPTOR-CHARGED_BETA")
    binana_features["catPi OTHER RECEPTOR"] = main_binana_out['pi_cation'].get("PI-CATION_RECEPTOR-CHARGED_OTHER")
    binana_features["catPi ALPHA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_ALPHA")
    binana_features["catPi BETA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binana_features["catPi OTHER LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")

    #5. ok add rotatable bond count to binana features
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
binana_df.to_csv('test_04_2.csv', index=False)
