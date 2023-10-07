# This script is designed to perform molecular docking 
# of decoys to their respective protein receptors.

from tqdm import tqdm
from joblib import Parallel, delayed
import os
import multiprocessing
from biopandas.pdb import PandasPdb
import pandas as pd
import shutil
import timeit

#  This script reads a CSV file containing molecular names and 
#  get the coordinates of the crystal ligand binding site.
def get_coordinates(file, padding): # file is the crystal ligand file
    # read the pdb file as a dataframe
    pmol = PandasPdb().read_pdb(file)

    # combine atoms and hetatm coordinates into one master dataframe
    atom_df = pmol.df['ATOM']
    hetatm_df = pmol.df['HETATM']
    full_atom_df = pd.concat([atom_df, hetatm_df])

    # find the minimum and maximum x, y, z coordinates for atoms in crystal ligand file
    x_min, x_max = min(full_atom_df.x_coord), max(full_atom_df.x_coord)
    y_min, y_max = min(full_atom_df.y_coord), max(full_atom_df.y_coord)
    z_min, z_max = min(full_atom_df.z_coord), max(full_atom_df.z_coord)

    # find center point for each axis
    x_center = x_min + abs(x_min-x_max)/2
    y_center = y_min + abs(y_min-y_max)/2
    z_center = z_min + abs(z_min-z_max)/2

    # calculate lengths for each axis based on difference between min and max values multiplied by user defined padding variable
    x_range = (abs(x_min-x_max)+padding)
    y_range = (abs(y_min-y_max)+padding)
    z_range = (abs(z_min-z_max)+padding)

    return x_center, y_center, z_center, x_range, y_range, z_range

# run the get_coordinates function in parallel and perform docking
# the docked files are saved in the docked_folder which are pdbqt format
def dock_file(docker_command, protein_filepath, ligand_filepath, center_x, center_y, center_z, size_x, size_y, size_z): # protein_filepath is the receptor file, ligand_filepath is the decoy file
    os.system(f'{docker_command} --receptor {protein_filepath} --ligand {ligand_filepath}  --center_x  {center_x} --center_y {center_y} --center_z {center_z} --size_x  {size_x} --size_y {size_y}  --size_z {size_z}' \
              ' --exhaustiveness=32 --num_wolves=40 --num_modes=5 --energy_range=4')
  
# read the CSV file containing molecular names and get the coordinates of the crystal ligand binding site
def dock_decoy(decoy_file, docker_command, receptor_folder, docked_folder, padding):
    filename = os.path.basename(decoy_file)
    pdb_code = filename.split('_')[0]
    receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
    example_crystal_ligand = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_ligand.pdbqt")

    try:
        # compute coordinates
        coordinates = get_coordinates(example_crystal_ligand, padding)
        print(coordinates)
        with open("filename.txt", "a") as bind_site_file:
            coordinates_str = ' '.join(map(str, coordinates))  # convert tuple elements to str and join with space
            bind_site_file.write(f"{filename} {coordinates_str}\n")  # use f-string

        # perform docking
        dock_file(docker_command, receptor_file, decoy_file, *coordinates)
        os.makedirs(docked_folder, exist_ok=True)
        origin_path = decoy_file.split('.')[0]
        shutil.move(f"{origin_path}_out.pdbqt", docked_folder)
        with open("filename_progess_rate.txt", "a") as progress_file:
            progress_file.write(f"{filename} docking finished\n")
        # shutil.copy(receptor_file, os.path.join(docked_folder, f"{pdb_code}_docked"))
        
    except Exception as e:
        with open("Errors.txt", "a") as error_file:  
            error_message = f"ERROR: Could not dock {filename}. Exception: {e}"
            print(error_message)
            error_file.write(error_message + "\n")  
    

# run the main function
start=timeit.default_timer()
docker_command = "filepath/SCORCH/utils/gwovina-1.0/build/linux/release/gwovina"
receptor_folder = "path to receptor folder"
decoy_folder = "path to decoy folder"
docked_folder = "path to docked folder"
padding = 12
decoy_files = [os.path.join(decoy_folder, file) for file in os.listdir(decoy_folder)]
with tqdm(total=len(decoy_files)) as pbar:
    Parallel(n_jobs=-1)(delayed(dock_decoy)(file, docker_command, receptor_folder, docked_folder, padding) for file in decoy_files)
    pbar.update()

end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))