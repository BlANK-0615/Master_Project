import random
import pandas as pd
import os
import shutil

def parse_pdbqt(file_path):
    """
    解析pdbqt文件，提取所有姿势。
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    poses = []
    pose = []
    for line in lines:
        if line.startswith("MODEL"):
            pose = []
        elif line.startswith("ENDMDL"):
            poses.append(pose)
            pose = []
        else:
            pose.append(line.strip())
            
    return poses

def random_select_and_record(poses, file_path, csv_path, output_folder):
    """
    对每个pdbqt文件随机选择一个姿势，并记录结果。
    """
    results = {}
    
    if len(poses) == 0:
        # 处理无效文件，记录错误信息
        error_message = "No valid poses found in the file"
        results[file_path] = error_message
        df = pd.DataFrame.from_dict(results, orient='index', columns=['Error']).T
        if not os.path.isfile(csv_path):
            df.to_csv(csv_path)
        else:
            df.to_csv(csv_path, mode='a', header=False)
        return
    
    selected_pose_index = random.randint(0, len(poses) - 1)
    
    for i in range(len(poses)):
        # 创建一个姿势文件的名称，格式为 pdbqt文件名+pose+数字
        pose_name = os.path.basename(file_path).replace("_out.pdbqt", "") + "_pose_" + str(i+1)

        # 判断是否为随机选择的姿势，如果是，则输出为独立的pdbqt文件，否则不输出
        if i == selected_pose_index:
            results[pose_name] = 1  # 选择的姿势
            output_path = os.path.join(output_folder, pose_name + ".pdbqt")
            with open(output_path, 'w') as output_file:
                output_file.write('\n'.join(poses[i]))
        else:
            results[pose_name] = 0  # 未选择的姿势
    
    # 将结果追加到CSV文件中
    df = pd.DataFrame.from_dict(results, orient='index', columns=['Selected']).T
    if not os.path.isfile(csv_path):
        df.to_csv(csv_path)
    else:
        df.to_csv(csv_path, mode='a', header=False)


def process_files(directory, output_folder, output_csv):
    """
    处理指定文件夹中的所有pdbqt文件，对每个文件进行随机选择姿势，并将结果记录到CSV文件中。
    同时将选中的姿势作为独立的pdbqt文件输出到另一个目录中。

    参数:
        directory: 包含pdbqt文件的文件夹的路径。
        output_folder: 输出独立姿势文件的目录路径。
        output_csv: 输出CSV文件的路径。
    """
    # 创建目标文件夹
    os.makedirs(output_folder, exist_ok=True)

    # 获取文件夹中的所有文件名
    file_names = os.listdir(directory)

    # 筛选出所有的pdbqt文件
    pdbqt_files = [f for f in file_names if f.endswith('.pdbqt')]

    for file_name in pdbqt_files:
        # 获取文件的完整路径
        file_path = os.path.join(directory, file_name)

        # 解析pdbqt文件，提取所有姿势
        poses = parse_pdbqt(file_path)

        # 对每个pdbqt文件随机选择一个姿势，并记录结果
        random_select_and_record(poses, file_path, output_csv, output_folder)


# 测试这个函数
process_files("/home/s2331261/Master_Project/6_Feature_Selection/test_get_pose_pdbqt", "/home/s2331261/Master_Project/6_Feature_Selection/tset_pose_results", "/home/s2331261/Master_Project/6_Feature_Selection/test_pose.csv")
