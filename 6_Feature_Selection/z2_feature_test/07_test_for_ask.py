# import os and set tensorflow verbosity
import os
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
# os.environ['NUMEXPR_MAX_THREADS'] = '1'

# import other libraries
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
    keep_closest_contacts = ['2.5 (HD, OA)', '2.5 (HD, HD)', '2.5 (C, HD)']

    #1.ok add closest contacts to binana_features dict(2.5)
    for contact in keep_closest_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ','_')
        binana_features[contact] = main_binana_out['closest'].get(binana_name)

    # return dictionary
    return binana_features


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


# 从pdbqt文件中读取配体姿态信息
ligand_file= "/home/s2331261/Master_Project/6_Feature_Selection/C1_get_decoy_poses/01_decoy_split_poses/1a0q_decoy_231_pose_3.pdbqt"
receptor_folder = "/home/s2331261/Master_Project/z1_3p_dataset/dock_source/all_receptors_ligands"
filename = os.path.basename(ligand_file)
pdb_code = filename.split('_')[0]
receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
pose_name = os.path.splitext(filename)[0]

lig_text = open(ligand_file, 'r').read()
lines = lig_text.split('\n')
clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
if len(clean_lines) < 3:
    pass
else:
    pose_info = '\n'.join(clean_lines)


# extract the interaction features
features = extract(pose_info,receptor_file)
print(features)

# fill None values with 0 for binana features
# multi_pose_features=features.fillna(0)
features.fillna(0, inplace=True)
multi_pose_features = features
print("multi_pose_features:\n",multi_pose_features)
multi_pose_features['Pose_name'] = pose_name
# binana_dict = run_binana(pose_info,receptor_file)
# binana_df = pd.DataFrame([binana_dict])
# print("#######################################################################")
# print(binana_df)
multi_pose_features.to_csv('test_06.csv', index=False)
