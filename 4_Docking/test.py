import os
import shutil

os.makedirs(docked_folder, exist_ok=True)
decoy_folder = "/home/s2331261/Master_Project/3_Docking/test_dock"
decoy_files = [os.path.join(decoy_folder, file) for file in os.listdir(decoy_folder)] 

filename = os.path.basename(decoy_file)
pdb_code = filename.split('_')[0]
foldername = filename.split('.')[0]
receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
example_crystal_ligand = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_ligand.pdbqt")

docked_folder = "eve"