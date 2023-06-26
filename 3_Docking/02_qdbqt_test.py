import os
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
from warnings import filterwarnings
import pprint
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from tqdm import tqdm
import shutil
import sys
from rdkit import RDLogger

# preparing a receptor
./utils/MGLTools-1.5.6/bin/pythonsh 
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py 
-r examples/predocked_1a0q/1a0q_receptor.pdb 
-A hydrogens \
-o examples/predocked_1a0q/1a0q_receptor_out.pdbqt 
-U nphs

def autodock_convert(destination_path, mgl_tools_path, ligand_filepath): # converts files from .pdb format to .pdbqt format using AutoDockTools

    # setup for populating
    errors_list = list()

    # define mgl_tools script paths
    pythonsh_path = f'{mgl_tools_path}bin/pythonsh'
    prepare_ligand_path = f'{mgl_tools_path}MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'

    # get filename from filepath
    ligand_file = ligand_filepath.split('/')[len(ligand_filepath.split('/')) - 1]

    # use AutoDockTools CLI to convert pdb file adding gasteiger charges and polar hydrogens
    try:
        output = subprocess.check_output(f'{pythonsh_path} {prepare_ligand_path} -l {ligand_filepath} -A hydrogens -o {destination_path}{ligand_file}qt -U nphs', shell=True, stderr=subprocess.STDOUT)
    except:
        print(f'ERROR: {ligand_filepath}')
# ./utils/MGLTools-1.5.6/bin/pythonsh \
# ./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
# -r /home/s2331261/Master_Project/3_Docking/m1_H.pdb \
# -A hydrogens \
# -o /home/s2331261/Master_Project/3_Docking/m1_H.pdbqt \
# -U nphs


./utils/MGLTools-1.5.6/bin/pythonsh \
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l /home/s2331261/Master_Project/3_Docking/m1_H.pdb \
-A hydrogens \
-o /home/s2331261/Master_Project/3_Docking/m1_H.pdbqt \
-U nphs


1tnh_decoy_323.pdb

./utils/MGLTools-1.5.6/bin/pythonsh \
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l /home/s2331261/Master_Project/3_Docking/pdb_files/1a4k_decoy_824.pdb \
-A hydrogens \
-o /home/s2331261/Master_Project/3_Docking/pdb_test_files/1a4k_decoy_824.pdbqt \
-U nphs

./utils/MGLTools-1.5.6/bin/pythonsh \
./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-l /home/s2331261/Master_Project/3_Docking/pdb_files/1tnh_decoy_323.pdb \
-A hydrogens \
-o /home/s2331261/Master_Project/3_Docking/pdb_test_files/1tnh_decoy_323.pdbqt \
-U nphs