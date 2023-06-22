import os
import shutil

# 指定要复制的文件数量
n = 10

# 创建目标文件夹

# 获取pdb_files文件夹中的文件列表
file_names = os.listdir('pdb_files')

# 复制前n个文件到pdb_test_files文件夹中
for file_name in file_names[:n]:
    src = os.path.join('pdb_files', file_name)
    dst = os.path.join('pdb_test_files', file_name)
    shutil.copy2(src, dst)
