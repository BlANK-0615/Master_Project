
# import os and set tensorflow verbosity
import os
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# os.environ['NUMEXPR_MAX_THREADS'] = '1'

# import other libraries
import sys
import glob # new module
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


# 1.Initialize logging
logging.basicConfig(filename='feature_log.txt', level=logging.INFO, format='%(asctime)s %(message)s')


# 2. Define binana features
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
    keep_closest_contacts = ["2.5 (HD, OA)","2.5 (HD, HD)","2.5 (HD, N)","2.5 (C, HD)","2.5 (MN, OA)","2.5 (HD, MN)","2.5 (OA, ZN)","2.5 (HD, ZN)","2.5 (A, HD)","2.5 (HD, NA)","2.5 (CA, HD)","2.5 (CA, OA)","2.5 (HD, MG)","2.5 (HD, P)","2.5 (F, HD)","2.5 (BR, HD)","2.5 (C, F)","2.5 (OA, OA)","2.5 (NA, OA)","2.5 (HD, SA)","2.5 (MG, OA)","2.5 (FE, HD)","2.5 (FE, OA)","2.5 (FE, N)","2.5 (A, N)","2.5 (N, N)","2.5 (N, OA)","2.5 (C, N)","2.5 (C, OA)","2.5 (CL, N)","2.5 (CL, HD)","2.5 (C, CL)","2.5 (A, A)","2.5 (A, C)","2.5 (A, OA)","2.5 (CL, OA)","2.5 (NA, ZN)","2.5 (C, C)","2.5 (C, SA)","2.5 (HD, S)","2.5 (HD, NI)","2.5 (NI, OA)","2.5 (CA, NA)","2.5 (C, MG)","2.5 (CO, OA)","2.5 (CO, HD)","2.5 (OA, SA)","2.5 (F, OA)","2.5 (C, NA)","2.5 (N, SA)","2.5 (F, ZN)","2.5 (F, N)","2.5 (A, F)","2.5 (MG, N)","2.5 (A, CL)","2.5 (C, ZN)","2.5 (N, ZN)","2.5 (MN, NA)","2.5 (MG, NA)","2.5 (A, NA)","2.5 (F, MG)","2.5 (K, OA)","2.5 (K, P)","2.5 (A, ZN)","2.5 (FE, NA)","2.5 (CU, OA)","2.5 (CU, NA)","2.5 (CU, HD)","2.5 (C, MN)","2.5 (N, NI)","2.5 (NA, NI)","2.5 (A, MG)","2.5 (CO, NA)","2.5 (MG, P)","2.5 (NA, NA)","2.5 (MN, N)","2.5 (CL, SA)","2.5 (A, CA)","2.5 (A, SA)","2.5 (OA, P)","2.5 (N, NA)","2.5 (F, SA)","2.5 (C, FE)","2.5 (A, FE)","2.5 (F, FE)","2.5 (MG, SA)","2.5 (N, P)","2.5 (C, CA)","2.5 (BR, N)","2.5 (BR, OA)","2.5 (C, CO)","2.5 (CO, N)","2.5 (A, MN)","2.5 (CL, FE)","2.5 (A, BR)","2.5 (SA, ZN)","2.5 (HD, K)","2.5 (CD, HD)","2.5 (CA, N)","2.5 (C, K)","2.5 (OA, S)","2.5 (HD, I)","2.5 (MN, P)","2.5 (C, P)","2.5 (FE, SA)","2.5 (F, MN)","2.5 (SA, SA)","2.5 (K, N)","2.5 (C, NI)","2.5 (BR, C)","2.5 (CD, OA)","2.5 (CL, MG)","2.5 (CL, ZN)","2.5 (S, ZN)","2.5 (NA, P)","2.5 (N, S)","2.5 (C, S)","2.5 (NA, SA)","2.5 (CU, N)","2.5 (A, P)","2.5 (A, CO)","2.5 (A, NI)","2.5 (C, CU)","2.5 (BR, MN)","2.5 (I, N)","2.5 (CL, MN)","2.5 (K, NA)","2.5 (BR, SA)","2.5 (A, K)","2.5 (CA, P)"]
    #len(keep_closest_contacts)
    # 130
    keep_close_contacts = ["4.0 (C, C)","4.0 (HD, OA)","4.0 (C, HD)","4.0 (C, N)","4.0 (A, C)","4.0 (A, OA)","4.0 (N, OA)","4.0 (A, N)","4.0 (HD, N)","4.0 (HD, HD)","4.0 (A, HD)","4.0 (OA, OA)","4.0 (C, OA)","4.0 (N, N)","4.0 (MN, OA)","4.0 (C, MN)","4.0 (C, SA)","4.0 (HD, SA)","4.0 (OA, SA)","4.0 (MN, SA)","4.0 (HD, MN)","4.0 (N, SA)","4.0 (A, A)","4.0 (N, S)","4.0 (HD, S)","4.0 (OA, S)","4.0 (S, ZN)","4.0 (N, ZN)","4.0 (HD, ZN)","4.0 (A, ZN)","4.0 (A, SA)","4.0 (SA, ZN)","4.0 (OA, ZN)","4.0 (C, S)","4.0 (A, S)","4.0 (C, ZN)","4.0 (C, NA)","4.0 (NA, OA)","4.0 (NA, SA)","4.0 (HD, NA)","4.0 (N, NA)","4.0 (A, NA)","4.0 (C, CA)","4.0 (CA, OA)","4.0 (CA, N)","4.0 (CA, HD)","4.0 (BR, C)","4.0 (BR, HD)","4.0 (BR, OA)","4.0 (BR, N)","4.0 (C, P)","4.0 (HD, P)","4.0 (OA, P)","4.0 (N, P)","4.0 (A, P)","4.0 (C, MG)","4.0 (MG, N)","4.0 (HD, MG)","4.0 (MG, OA)","4.0 (MG, P)","4.0 (C, F)","4.0 (F, N)","4.0 (F, OA)","4.0 (F, HD)","4.0 (A, F)","4.0 (NA, ZN)","4.0 (CL, OA)","4.0 (C, CL)","4.0 (CL, N)","4.0 (CL, HD)","4.0 (A, CL)","4.0 (A, BR)","4.0 (MN, N)","4.0 (A, MN)","4.0 (F, MN)","4.0 (CL, SA)","4.0 (BR, SA)","4.0 (CA, SA)","4.0 (C, FE)","4.0 (FE, OA)","4.0 (FE, N)","4.0 (FE, HD)","4.0 (F, SA)","4.0 (A, CA)","4.0 (SA, SA)","4.0 (C, NI)","4.0 (MN, P)","4.0 (F, ZN)","4.0 (HD, NI)","4.0 (NI, S)","4.0 (NI, OA)","4.0 (N, NI)","4.0 (MG, NA)","4.0 (CA, P)","4.0 (CA, NA)","4.0 (A, MG)","4.0 (CL, ZN)","4.0 (P, ZN)","4.0 (CO, P)","4.0 (C, CO)","4.0 (CO, HD)","4.0 (CO, OA)","4.0 (CA, S)","4.0 (NA, NI)","4.0 (P, SA)","4.0 (CL, MG)","4.0 (CO, SA)","4.0 (FE, NA)","4.0 (FE, P)","4.0 (A, CU)","4.0 (C, CU)","4.0 (A, FE)","4.0 (MN, NA)","4.0 (CL, FE)","4.0 (F, FE)","4.0 (CU, NA)","4.0 (K, P)","4.0 (K, OA)","4.0 (C, K)","4.0 (HD, K)","4.0 (MG, SA)","4.0 (S, SA)","4.0 (FE, SA)","4.0 (A, CO)","4.0 (CO, N)","4.0 (CU, N)","4.0 (CU, HD)","4.0 (CU, OA)","4.0 (MG, S)","4.0 (NI, SA)","4.0 (A, NI)","4.0 (CL, NI)","4.0 (BR, MG)","4.0 (CO, NA)","4.0 (BR, MN)","4.0 (NA, NA)","4.0 (HD, I)","4.0 (C, I)","4.0 (I, OA)","4.0 (BR, ZN)","4.0 (F, MG)","4.0 (NA, P)","4.0 (CL, MN)","4.0 (CO, S)","4.0 (CA, F)","4.0 (C, CD)","4.0 (I, N)","4.0 (CU, F)","4.0 (FE, S)","4.0 (A, I)","4.0 (I, ZN)","4.0 (I, SA)","4.0 (A, K)","4.0 (K, N)","4.0 (CD, P)","4.0 (CD, OA)","4.0 (CD, HD)","4.0 (MN, S)","4.0 (F, K)","4.0 (CL, NA)","4.0 (CL, CU)","4.0 (CO, F)","4.0 (K, SA)","4.0 (K, NA)","4.0 (CU, SA)","4.0 (A, CD)","4.0 (CD, N)","4.0 (CL, CO)","4.0 (CA, CL)","4.0 (F, NA)","4.0 (HG, OA)","4.0 (HD, HG)","4.0 (NA, S)","4.0 (I, MG)"]
    #len(keep_close_contacts)
    # 174
    keep_ligand_atoms = ['LA C','LA OA','LA N','LA HD','LA SA','LA A','LA S','LA NA','LA BR','LA P','LA F','LA CL','LA I']
    #len(keep_ligand_atoms)
    #13
    keep_elsums = ["ElSum (C, C)","ElSum (HD, OA)","ElSum (C, HD)","ElSum (C, N)","ElSum (A, C)","ElSum (A, OA)","ElSum (N, OA)","ElSum (A, N)","ElSum (HD, N)","ElSum (HD, HD)","ElSum (A, HD)","ElSum (OA, OA)","ElSum (C, OA)","ElSum (N, N)","ElSum (MN, OA)","ElSum (C, MN)","ElSum (C, SA)","ElSum (HD, SA)","ElSum (OA, SA)","ElSum (MN, SA)","ElSum (HD, MN)","ElSum (N, SA)","ElSum (A, A)","ElSum (N, S)","ElSum (HD, S)","ElSum (OA, S)","ElSum (S, ZN)","ElSum (OA, ZN)","ElSum (N, ZN)","ElSum (HD, ZN)","ElSum (A, ZN)","ElSum (A, SA)","ElSum (SA, ZN)","ElSum (C, S)","ElSum (A, S)","ElSum (C, ZN)","ElSum (C, NA)","ElSum (NA, OA)","ElSum (NA, SA)","ElSum (HD, NA)","ElSum (N, NA)","ElSum (A, NA)","ElSum (C, CA)","ElSum (CA, OA)","ElSum (CA, N)","ElSum (CA, HD)","ElSum (BR, C)","ElSum (BR, HD)","ElSum (BR, OA)","ElSum (BR, N)","ElSum (C, P)","ElSum (HD, P)","ElSum (OA, P)","ElSum (N, P)","ElSum (A, P)","ElSum (C, MG)","ElSum (MG, N)","ElSum (HD, MG)","ElSum (MG, OA)","ElSum (MG, P)","ElSum (C, F)","ElSum (F, N)","ElSum (F, HD)","ElSum (F, OA)","ElSum (A, F)","ElSum (NA, ZN)","ElSum (CL, OA)","ElSum (C, CL)","ElSum (CL, N)","ElSum (CL, HD)","ElSum (A, CL)","ElSum (A, BR)","ElSum (MN, N)","ElSum (A, MN)","ElSum (F, MN)","ElSum (CL, SA)","ElSum (BR, SA)","ElSum (CA, SA)","ElSum (C, FE)","ElSum (FE, OA)","ElSum (FE, N)","ElSum (FE, HD)","ElSum (F, SA)","ElSum (A, CA)","ElSum (SA, SA)","ElSum (C, NI)","ElSum (MN, P)","ElSum (F, ZN)","ElSum (HD, NI)","ElSum (NI, S)","ElSum (NI, OA)","ElSum (N, NI)","ElSum (MG, NA)","ElSum (CA, P)","ElSum (CA, NA)","ElSum (A, MG)","ElSum (CL, ZN)","ElSum (P, ZN)","ElSum (CO, P)","ElSum (C, CO)","ElSum (CO, OA)","ElSum (CO, HD)","ElSum (CA, S)","ElSum (NA, NI)","ElSum (P, SA)","ElSum (CL, MG)","ElSum (CO, SA)","ElSum (FE, NA)","ElSum (FE, P)","ElSum (A, CU)","ElSum (C, CU)","ElSum (A, FE)","ElSum (MN, NA)","ElSum (CL, FE)","ElSum (F, FE)","ElSum (CU, NA)","ElSum (F, MG)","ElSum (K, P)","ElSum (K, OA)","ElSum (C, K)","ElSum (HD, K)","ElSum (MG, SA)","ElSum (S, SA)","ElSum (FE, SA)","ElSum (A, CO)","ElSum (CO, N)","ElSum (CU, N)","ElSum (CU, HD)","ElSum (CU, OA)","ElSum (MG, S)","ElSum (NI, SA)","ElSum (A, NI)","ElSum (CL, NI)","ElSum (BR, MG)","ElSum (CO, NA)","ElSum (BR, MN)","ElSum (NA, NA)","ElSum (HD, I)","ElSum (C, I)","ElSum (I, OA)","ElSum (BR, ZN)","ElSum (NA, P)","ElSum (CL, MN)","ElSum (CO, S)","ElSum (CA, F)","ElSum (C, CD)","ElSum (I, N)","ElSum (CU, F)","ElSum (FE, S)","ElSum (A, I)","ElSum (I, ZN)","ElSum (I, SA)","ElSum (A, K)","ElSum (K, N)","ElSum (CD, P)","ElSum (CD, OA)","ElSum (CD, HD)","ElSum (MN, S)","ElSum (F, K)","ElSum (CL, NA)","ElSum (CL, CU)","ElSum (CO, F)","ElSum (K, SA)","ElSum (K, NA)","ElSum (CU, SA)","ElSum (A, CD)","ElSum (CD, N)","ElSum (CL, CO)","ElSum (CA, CL)","ElSum (F, NA)","ElSum (HG, OA)","ElSum (HD, HG)","ElSum (NA, S)","ElSum (I, MG)"]
    #len(keep_elsums)
    #174
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
    binana_features["catPi ALPHA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_ALPHA")
    binana_features["catPi BETA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binana_features["catPi OTHER LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")
    binana_features["catPi ALPHA RECEPTOR"] = main_binana_out['pi_cation'].get("PI-CATION_RECEPTOR-CHARGED_ALPHA")
    binana_features["catPi BETA RECEPTOR"] = main_binana_out['pi_cation'].get("PI-CATION_RECEPTOR-CHARGED_BETA")
    binana_features["catPi OTHER RECEPTOR"] = main_binana_out['pi_cation'].get("PI-CATION_RECEPTOR-CHARGED_OTHER")

    #5. ok add rotatable bond count to binana features
    binana_features["nRot"] = main_binana_out['nrot']

    # return dictionary
    return binana_features

# 3. Define kier flexibility
def kier_flexibility(pose_info):

    """
    Function: Calculate Kier flexibility    
    for ligand                              
                                            
    Inputs: ligand as a pdbqt string block  
                                            
    Output: Kier flexibility                
    """

    # parse pdbqt block
    mol = kier.SmilePrep(pose_info)

    # calculate flexibility
    return kier.CalculateFlexibility(mol)

# 4. Define ECIFs
def calculate_ecifs(pose_info,receptor_file):

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
    ECIF_data = GetECIF(receptor_file, pose_info, distance_cutoff=6.0)

    # replace the semicolons to make valid dataframe headers
    ECIFHeaders = [header.replace(';','') for header in PossibleECIF]

    # zip into a dictionary and convert to dataframe
    ECIF_data = dict(zip(ECIFHeaders,ECIF_data))
    ECIF_df = pd.DataFrame(ECIF_data,index=[0])

    # return the dataframe
    return ECIF_df


#5. Define extract function
def extract(pose_info,receptor_file):

    """
    Function: Get all descriptor features   
    for protein-ligand complex              
                                            
    Inputs: ligand as a pdbqt string block, 
    receptor pdbqt filepath     
                                            
    Output: All protein-ligand complex      
    descriptor features as a DataFrame      
    """
    
    # get the kier flexibility
    k = kier_flexibility(pose_info)

    # get the binana descriptors and build into dataframe
    binana_dict = run_binana(pose_info,receptor_file)
    binana_df = pd.DataFrame([binana_dict])

    # get the ECIFs
    ECIF = calculate_ecifs(pose_info,receptor_file)

    # concatenate all feature columns to one row dataframe
    df = pd.concat([ECIF,binana_df],axis=1)

    # add the kier flexibility as a column
    df['Kier Flexibility'] = k

    # return the features
    return df


# 6. Define extract_features function
def extract_features(ligand_file, receptor_file):
    # your original code...
    filename = os.path.basename(ligand_file)
    pdb_code = filename.split('_')[0]
    pose_name = os.path.splitext(filename)[0]

    lig_text = open(ligand_file, 'r').read()
    lines = lig_text.split('\n')
    clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
    if len(clean_lines) < 3:
        logging.warning(f"Failed to process ligand file {ligand_file}")
        return None
    else:
        pose_info = '\n'.join(clean_lines)
    
    # extract the interaction features
    features = extract(pose_info, receptor_file)
    features.fillna(0, inplace=True)
    multi_pose_features = features
    multi_pose_features['Pose_name'] = pose_name

    return multi_pose_features


def process_all(receptor_folder, split_ligand_folder):
    # Find all subfolders
    # folders = [f.path for f in os.scandir(protein_ligand_directory) if f.is_dir()]
    
    # Create an empty DataFrame to store all features
    all_features = pd.DataFrame()
    # get the list of ligand files
    split_files = [os.path.join(split_ligand_folder, file) for file in os.listdir(split_ligand_folder)]
    
    # run the extract_features function on each ligand file
    for split_file in split_files:
        # get the filename
        filename = os.path.basename(split_file)
        pdb_code = filename.split('_')[0]
        # get the receptor file and ligand file
        receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
        ligand_file = split_file
        features = extract_features(ligand_file, receptor_file) 
        if features is not None:
            all_features = all_features.append(features)
            logging.info(f"Processing ligand: {ligand_file} and receptor: {receptor_file}")

    # for folder in folders:
    #     ligand_file = glob.glob(os.path.join(folder, '*_ligand.pdbqt'))[0]
    #     receptor_file = glob.glob(os.path.join(folder, '*_receptor.pdbqt'))[0]       
    #     features = extract_features(ligand_file, receptor_file)      
    #     if features is not None:
    #         all_features = all_features.append(features)
    #         logging.info(f"Processing ligand: {ligand_file} and receptor: {receptor_file}")
    
    # Save all features to a single csv file
    all_features.to_csv('all_features_0716_01.csv', index=False)
    # Save all features to a hdf5 file
    # all_features.to_hdf('all_features.h5', key='df', mode='w')

# MAIN 
# define the receptor and ligand folders
receptor_folder = "/home/s2331261/Master_Project/z1_3p_dataset/dock_source/all_receptors_ligands"
split_ligand_folder = "/home/s2331261/Master_Project/6_Feature_Selection/C1_get_decoy_poses/01_decoy_split_poses"
# Call the function with the directory containing all the proteins and ligands
process_all(receptor_folder, split_ligand_folder)


# # get the list of ligand files
# split_files = [os.path.join(split_ligand_folder, file) for file in os.listdir(split_ligand_folder)]
# for split_file in split_files:
#     # get the filename
#     filename = os.path.basename(split_file)
#     pdb_code = filename.split('_')[0]
#     # get the receptor file
#     receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
#     ligand_file = split_file


# def dock_decoy(decoy_file, docker_command, receptor_folder, docked_folder, padding):
#     filename = os.path.basename(decoy_file)
#     pdb_code = filename.split('_')[0]
#     receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
#     example_crystal_ligand = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_ligand.pdbqt")