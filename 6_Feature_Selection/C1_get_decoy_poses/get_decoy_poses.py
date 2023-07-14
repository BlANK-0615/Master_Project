import random
import pandas as pd
import os
import shutil
import timeit


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
            pose = [line.strip()]  # 添加MODEL行到姿势中
        elif line.startswith("ENDMDL"):
            pose.append(line.strip())  # 添加ENDMDL行到姿势中
            poses.append(pose)
            pose = []
        else:
            pose.append(line.strip())
            
    return poses

def random_select_and_record(poses, file_path, csv_path, output_folder):
    """
    对每个pdbqt文件随机选择一个姿势，并记录结果到CSV文件。
    """
    results = {}
    
    if len(poses) == 0:
        # 处理无效文件，记录错误信息
        error_message = "No valid poses found in the file"
        results[file_path] = error_message
        with open("Error_poses.txt", 'a') as error_file:
            error_file.write(f"{file_path}: {error_message}\n")
    
    selected_pose_index = random.randint(0, len(poses) - 1)
    
    selected_pose = poses[selected_pose_index]
    
    pose_name = os.path.basename(file_path).replace("_out.pdbqt", "")
    pose_no = selected_pose_index + 1
    
    output_path = os.path.join(output_folder, f"{pose_name}_pose_{pose_no}.pdbqt")
    with open(output_path, 'w') as output_file:
        output_file.write('\n'.join(selected_pose))

    results["Pose name"] = pose_name
    results["Pose_no"] = pose_no
    results["Label"] = 0
    
    df = pd.DataFrame(results, index=[0])
    if not os.path.isfile(csv_path):
        df.to_csv(csv_path, index=False)
    else:
        df.to_csv(csv_path, mode='a', header=False, index=False)


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

        # 对每个pdbqt文件随机选择一个姿势，并记录结果到CSV文件
        random_select_and_record(poses, file_path, output_csv, output_folder)


# 测试这个函数
start=timeit.default_timer()
process_files("/home/s2331261/Master_Project/6_Feature_Selection/A3_nop_docked_pdbqt",
                "/home/s2331261/Master_Project/6_Feature_Selection/C1_get_decoy_poses/01_decoy_split_poses", 
                "/home/s2331261/Master_Project/6_Feature_Selection/C1_get_decoy_poses/decoy_log.csv")
# process_files("源文件夹路径", "输出文件夹路径", "输出CSV文件路径")
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))
