import os
import shutil

# 源文件夹和目标文件夹
src_folder = '/home/s2331261/Master_Project/3_Docking/A2_pdbqt_files'
dst_folder = '/home/s2331261/Master_Project/3_Docking/A4_deep_pdbqt_files'

# 确保目标文件夹存在
os.makedirs(dst_folder, exist_ok=True)

# 获取源文件夹中的所有文件
file_names = os.listdir(src_folder)

# 遍历文件
for file_name in file_names:
    # 检查文件名是否符合指定的模式
    if '_deep' in file_name:
        # 计算源文件和目标文件的全路径
        src = os.path.join(src_folder, file_name)
        dst = os.path.join(dst_folder, file_name)
        
        # 移动文件
        shutil.move(src, dst)
        # print(src, '->', dst)
