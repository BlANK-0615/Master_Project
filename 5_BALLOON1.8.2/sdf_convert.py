import pandas as pd
import os
import subprocess
from joblib import Parallel, delayed
from tqdm import tqdm

# 错误日志文件
sdf_error_file = "/home/s2331261/Master_Project/BALLOON1.8.2/sdf_error_logs.txt"

def generate_sdf(smiles, name, dir_path):
    try:
        # 设置文件名
        name_new = name.replace("_", "_decoy_") + "_rdm"
        sdf_file = os.path.join(dir_path, f'{name_new}.sdf')

        # 使用balloon将SMILES转换为SDF
        balloon_cmd = [
            '/home/s2331261/Master_Project/BALLOON1.8.2/balloon',
            '-f', '/home/s2331261/Master_Project/BALLOON1.8.2/MMFF94.mff',
            '--nconfs', '1', 
            '--noGA', smiles, sdf_file
        ]
        subprocess.run(balloon_cmd, check=True)
    except Exception as e:
        with open(sdf_error_file, "a") as err:
            err.write(f"Error processing {name}: {str(e)}\n")

# 读取CSV文件
df = pd.read_csv('/home/s2331261/Master_Project/3_Docking/B2_results_files/y1_broken_deocy_emb(total).csv')

# 为每个SMILES生成SDF文件
Parallel(n_jobs=-1)(delayed(generate_sdf)(row['smile'], row['name'], '/home/s2331261/Master_Project/BALLOON1.8.2/A1_correct_sdf_files') for i, row in tqdm(df.iterrows(), total=df.shape[0]))


# /home/s2331261/Master_Project/BALLOON1.8.2/balloon -f /home/s2331261/Master_Project/BALLOON1.8.2/MMFF94.mff \
#     --nconfs 1 --noGA "CC1CCC2C(C)(C)C(O)(C(=O)N=CC=CC(=O)CC3=NC(N[N+](C)=O)=NC3)C3(C)CCC14CC21CCC(O)C1(C)CC43" 5u6g.sdf