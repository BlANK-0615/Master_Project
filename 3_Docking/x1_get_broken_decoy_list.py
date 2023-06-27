import pandas as pd
import os

#df = pd.read_csv('/home/s2331261/Master_Project/3_Docking/Broken_Decoy.txt', sep=" ")

# 读取 TXT 文件，不读入列名
df = pd.read_csv('/home/s2331261/Master_Project/3_Docking/a8_broken.txt', sep=" ", header=None)

# 去除第五列中的冒号
df[4] = df[4].str.rstrip(':')

# 选择第五列和第六列，存储到新的 DataFrame
new_df = df[[4, 5]].copy()

# 添加列名
new_df.columns = ['name', 'smile']

# 保存到 CSV 文件
new_df.to_csv('a8_broken_deocy.csv', index=False)



