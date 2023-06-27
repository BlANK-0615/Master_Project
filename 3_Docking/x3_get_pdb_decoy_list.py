import pandas as pd
import os

folder_path = '/home/s2331261/Master_Project/3_Docking'  # 替换为你的文件夹路径

# 获取文件夹中的所有文件名
file_names = os.listdir(folder_path)

# 去除文件名中的 "decoy" 并保留格式为 "1a0q_231" 的部分
#file_names = [name.replace('_decoy', '').split('.')[0] for name in file_names]

# 创建 DataFrame 并将文件名存储到一列中
df = pd.DataFrame({'name': file_names})

# 保存到 CSV 文件
df.to_csv('3_Docking_list.csv', index=False)

# ls -l | grep "^-" | wc -l
