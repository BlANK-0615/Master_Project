import os
from rdkit import Chem
from joblib import Parallel, delayed
from tqdm import tqdm

# 错误日志文件
pdb_error_file = "/home/s2331261/Master_Project/BALLOON1.8.2/pdb_error_logs.txt"

# SDF和PDB文件路径
dir_path_sdf = '/home/s2331261/Master_Project/BALLOON1.8.2/A1_correct_sdf_files'
dir_path_pdb = '/home/s2331261/Master_Project/BALLOON1.8.2/A2_correct_pdb_files'

def convert_sdf_to_pdb(sdf_filename):
    try:
        name = os.path.splitext(sdf_filename)[0]
        sdf_file = os.path.join(dir_path_sdf, sdf_filename)
        pdb_file = os.path.join(dir_path_pdb, f'{name}.pdb')

        # 读取SDF文件
        mol_supplier = Chem.SDMolSupplier(sdf_file)
        mol = mol_supplier[0]

        # 移除氢原子并保存到PDB文件
        mol_noH = Chem.RemoveHs(mol)
        Chem.MolToPDBFile(mol_noH, pdb_file)
    except Exception as e:
        with open(pdb_error_file, "a") as err:
            err.write(f"Error processing {sdf_filename}: {str(e)}\n")

# 获取所有SDF文件名
sdf_files = [f for f in os.listdir(dir_path_sdf) if f.endswith('.sdf')]
# 并行转换所有SDF到PDB
Parallel(n_jobs=-1)(delayed(convert_sdf_to_pdb)(sdf_file) for sdf_file in tqdm(sdf_files, total=len(sdf_files)))