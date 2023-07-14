# import os
# import timeit


# def count_conformations(file):
#     with open(file, 'r') as f:
#         lines = f.readlines()
#     count = sum(1 for line in lines if line.find('MODEL') != -1)
#     return count

# def main():
#     folder_path = "/home/s2331261/Master_Project/4_Docking/A2_random_docked/random_decoy_docked" # 你的文件夹路径
#     filenames = os.listdir(folder_path) # 得到文件夹中所有文件的名称
#     pdbqt_files = [file for file in filenames if file.endswith('.pdbqt')] # 选取所有的pdbqt文件
    
#     results = {}
#     for file in pdbqt_files:
#         full_path = os.path.join(folder_path, file) # 获取文件的完整路径
#         num_conformations = count_conformations(full_path) # 计算每个文件中的构象数量
#         results[file] = num_conformations

#     # 将结果写入文件
#     with open('conformations_count.txt', 'w') as f:
#         for file, count in results.items():
#             f.write(f"{file}: {count}\n")

# if __name__ == "__main__":
#     start=timeit.default_timer()
#     main()
#     end=timeit.default_timer()
#     print('Running time: %s Seconds'%(end-start))

import os
import timeit

def count_conformations(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    count = sum(1 for line in lines if line.find('MODEL') != -1)
    return count

def main():
    folder_path = "/home/s2331261/Master_Project/4_Docking/A1_deep_docked/deep_decoy_docked" # 你的文件夹路径
    filenames = os.listdir(folder_path) # 得到文件夹中所有文件的名称
    pdbqt_files = [file for file in filenames if file.endswith('.pdbqt')] # 选取所有的pdbqt文件

    results = {}
    exception_results = {} # 创建一个新的字典来保存异常的结果
    for file in pdbqt_files:
        full_path = os.path.join(folder_path, file) # 获取文件的完整路径
        num_conformations = count_conformations(full_path) # 计算每个文件中的构象数量
        results[file] = num_conformations
        if num_conformations != 5: # 检查构象的数量是否等于5
            exception_results[file] = num_conformations

    # 将结果写入文件
    with open('CC.txt', 'a') as f:
        for file, count in results.items():
            f.write(f"{file}: {count}\n")

    # 将异常结果写入另一个文件
    with open('Error_CC.txt', 'a') as f:
        for file, count in exception_results.items():
            f.write(f"{file}: {count}\n")

if __name__ == "__main__":
    start=timeit.default_timer()
    main()
    end=timeit.default_timer()
    print('Running time: %s Seconds'%(end-start))








