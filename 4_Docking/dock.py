from tqdm import tqdm
from joblib import Parallel, delayed
import os
import multiprocessing
from biopandas.pdb import PandasPdb
import pandas as pd
import shutil
import timeit

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


def dock_file(docker_command, protein_filepath, ligand_filepath, center_x, center_y, center_z, size_x, size_y, size_z): # protein_filepath is the receptor file, ligand_filepath is the decoy file
    os.system(f'{docker_command} --receptor {protein_filepath} --ligand {ligand_filepath}  --center_x  {center_x} --center_y {center_y} --center_z {center_z} --size_x  {size_x} --size_y {size_y}  --size_z {size_z}' \
              ' --exhaustiveness=32 --num_wolves=40 --num_modes=5 --energy_range=4')
  
    
def dock_decoy(decoy_file, docker_command, receptor_folder, docked_folder, padding):
    filename = os.path.basename(decoy_file)
    pdb_code = filename.split('_')[0]
    receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")
    example_crystal_ligand = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_ligand.pdbqt")

    try:
        # compute coordinates
        coordinates = get_coordinates(example_crystal_ligand, padding)
        print(coordinates)
        with open("Toco_02_bind_site.txt", "a") as bind_site_file:
            coordinates_str = ' '.join(map(str, coordinates))  # convert tuple elements to str and join with space
            bind_site_file.write(f"{filename} {coordinates_str}\n")  # use f-string

        # perform docking
        dock_file(docker_command, receptor_file, decoy_file, *coordinates)
        os.makedirs(docked_folder, exist_ok=True)
        origin_path = decoy_file.split('.')[0]
        shutil.move(f"{origin_path}_out.pdbqt", docked_folder)
        with open("Toco_02_progess_rate.txt", "a") as progress_file:
            progress_file.write(f"{filename} docking finished\n")
        # shutil.copy(receptor_file, os.path.join(docked_folder, f"{pdb_code}_docked"))
        
    except Exception as e:
        with open("Errors.txt", "a") as error_file:  # "a" 表示追加模式，新的内容会被写入到文件的末尾
            error_message = f"ERROR: Could not dock {filename}. Exception: {e}"
            print(error_message)
            error_file.write(error_message + "\n")  # 添加一个换行符以便在文件中每条错误信息占一行
    

# def main(docker_command, receptor_folder, decoy_folder, docked_folder, padding):
#     decoy_files = [os.path.join(decoy_folder, file) for file in os.listdir(decoy_folder)] #file is the deocy file name
#     # generate a list of all decoy files
#     # num_cores = multiprocessing.cpu_count()
    
#     # with tqdm(total=len(decoy_files)) as pbar:
#     # def update(*a):
#     #     pbar.update()
#     Parallel(n_jobs=32)(delayed(dock_decoy)(file, docker_command, receptor_folder, docked_folder, padding) for file in tqdm(file, total=len(decoy_files)))
# #pbar.update()

#Parallel(n_jobs=32)(delayed(convert_pdb_to_pdbqt)(pdb_file) for pdb_file in tqdm(pdb_files, total=len(pdb_files)))



# run the main function
start=timeit.default_timer()

docker_command = "/home/s2331261/Master_Project/3_Docking/SCORCH/utils/gwovina-1.0/build/linux/release/gwovina"
receptor_folder = "/home/s2331261/Master_Project/z1_3p_dataset/dock_source/all_receptors_ligands"
decoy_folder = "/home/s2331261/Master_Project/3_Docking/decoy_qdbqt_02"
docked_folder = "Toco_docked_02"
padding = 12
decoy_files = [os.path.join(decoy_folder, file) for file in os.listdir(decoy_folder)]
# generate a list of all decoy files
# num_cores = multiprocessing.cpu_count()
with tqdm(total=len(decoy_files)) as pbar:
    # def update(*a):
    #     pbar.update()
    Parallel(n_jobs=-1)(delayed(dock_decoy)(file, docker_command, receptor_folder, docked_folder, padding) for file in decoy_files)
    pbar.update()
#Parallel(n_jobs=32)(delayed(convert_pdb_to_pdbqt)(pdb_file) for pdb_file in tqdm(pdb_files, total=len(pdb_files)))

# Parallel(n_jobs=-1)(delayed(dock_decoy)(file, docker_command, receptor_folder, docked_folder, padding) for file in tqdm(decoy_files, total=len(decoy_files)))

end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))