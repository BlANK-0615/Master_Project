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
        print(receptor_file)
        print(ligand_file)
        # features = extract_features(ligand_file, receptor_file) 
        # if features is not None:
        #     all_features = all_features.append(features)
        #     logging.info(f"Processing ligand: {ligand_file} and receptor: {receptor_file}")

receptor_folder = "/home/s2331261/Master_Project/z1_3p_dataset/dock_source/all_receptors_ligands"
split_ligand_folder = "/home/s2331261/Master_Project/6_Feature_Selection/C1_get_decoy_poses/02_decoy_test/test_pose_results_05"
# Call the function with the directory containing all the proteins and ligands
process_all(receptor_folder, split_ligand_folder)



