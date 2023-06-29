import os
import shutil

# 源文件夹和文件列表
src_folder = '/home/s2331261/Master_Project/3_Docking/A2_pdbqt_files'
file_names = os.listdir(src_folder)

# 文件夹名字前缀
folder_prefix = '/home/s2331261/Master_Project/3_Docking/decoy_qdbqt_'

# 将文件夹中的文件分成五份
chunks = [file_names[x:x+10000] for x in range(0, len(file_names), 10000)]

# 复制文件到各自的文件夹
for i, chunk in enumerate(chunks):
    # 创建目标文件夹
    dst_folder = folder_prefix + str(i+1).zfill(2)  # 在文件夹名字后添加数字，用0填充到两位数
    os.makedirs(dst_folder, exist_ok=True)

    # 复制文件
    for file_name in chunk:
        src = os.path.join(src_folder, file_name)
        dst = os.path.join(dst_folder, file_name)
        shutil.copy2(src, dst)


