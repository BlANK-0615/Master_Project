import os
import subprocess
from tqdm import tqdm

# 创建目标目录
# os.makedirs('pdbqt_files', exist_ok=True)

# # 列出pdb文件
# pdb_files = [f for f in os.listdir('pdb_files') if f.endswith('.pdb')]

# # 遍历所有pdb文件
# for i, pdb_file in enumerate(tqdm(pdb_files, total=len(pdb_files))):
#     # 获取文件的基本名（无后缀）
#     basename = os.path.splitext(pdb_file)[0]
#     base_parts = basename.split('_')
#     new_base = f'{base_parts[0]}_decoy_{base_parts[1]}'

#     # 构造输入文件和输出文件的路径
#     input_file = os.path.join('pdb_files', pdb_file)
#     output_file = os.path.join('pdbqt_files', f'{new_base}.pdbqt')

#     # 运行转换命令
#     command = [
#         './utils/MGLTools-1.5.6/bin/pythonsh',
#         './utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py',
#         '-l', input_file,
#         '-A', 'hydrogens',
#         '-o', output_file,
#         '-U', 'nphs',
#     ]
#     subprocess.run(command, check=True)


# 创建目标目录
os.makedirs('pdbqt_test_files', exist_ok=True)

# 列出pdb文件
pdb_files = os.listdir('pdb_test_files')

# 遍历所有pdb文件

for i, pdb_file in enumerate(tqdm(pdb_files, total=len(pdb_files))):
    if pdb_file.endswith('.pdb'):
        # 获取文件的基本名（无后缀）
        basename = os.path.splitext(pdb_file)[0]
        # 构造输入文件和输出文件的路径
        input_file = os.path.join('pdb_test_files', pdb_file)
        output_file = os.path.join('pdbqt_test_files', f'{basename}.pdbqt')
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
            with open('Broken_PDB.txt', 'a') as f:
                f.write(f"Error pdb file at index {pdb_file}\n")