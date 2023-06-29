import os
import subprocess
from joblib import Parallel, delayed
from tqdm import tqdm
import timeit

def convert_pdb_to_pdbqt(pdb_file):
    if pdb_file.endswith('.pdb'):
        # 获取文件的基本名（无后缀）
        basename = os.path.splitext(pdb_file)[0]
        # 构造输入文件和输出文件的路径
        input_file = os.path.join('/home/s2331261/Master_Project/3_Docking/end_broken_file/random_pdb_files', pdb_file)
        output_file = os.path.join('/home/s2331261/Master_Project/3_Docking/end_broken_file/random_pdbqt_files', f'{basename}_rdm.pdbqt')
        try:
            # 运行转换命令
            command = [
                '/home/s2331261/Master_Project/3_Docking/SCORCH/utils/MGLTools-1.5.6/bin/pythonsh',
                '/home/s2331261/Master_Project/3_Docking/SCORCH/utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py',
                '-l', input_file,
                '-A', 'hydrogens',
                '-o', output_file,
                '-U', 'nphs',
            ]
            subprocess.run(command, check=True)
        except:
            print(f'ERROR: {pdb_file}')
            with open('/home/s2331261/Master_Project/3_Docking/end_broken_file/pdbqt_broken.txt', 'a') as f:
                f.write(f"Error pdb file at index {pdb_file}\n")

# 创建目标目录
os.makedirs('/home/s2331261/Master_Project/3_Docking/end_broken_file/random_pdbqt_files', exist_ok=True)

# 列出pdb文件
pdb_files = os.listdir('/home/s2331261/Master_Project/3_Docking/end_broken_file/random_pdb_files')

# 并行转换所有pdb文件
# n_jobs参数可以根据你的实际情况进行调整，例如你的物理核心数或逻辑核心数
start=timeit.default_timer()
Parallel(n_jobs=32)(delayed(convert_pdb_to_pdbqt)(pdb_file) for pdb_file in tqdm(pdb_files, total=len(pdb_files)))
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))