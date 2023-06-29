import os
import shutil

# 指定要复制的文件数量
#n = 32

# 创建目标文件夹

# 获取pdb_files文件夹中的文件列表
file_names = os.listdir('/home/s2331261/Master_Project/3_Docking/A2_pdbqt_files')
os.makedirs('/home/s2331261/Master_Project/3_Docking/final_test_02', exist_ok=True)

# 复制前n个文件到pdb_test_files文件夹中
for file_name in file_names[32:64]:
    src = os.path.join('/home/s2331261/Master_Project/3_Docking/A2_pdbqt_files', file_name)
    dst = os.path.join('/home/s2331261/Master_Project/3_Docking/final_test_02', file_name)
    shutil.copy2(src, dst)
