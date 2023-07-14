import csv
import os
import shutil
import timeit


start=timeit.default_timer()
# CSV文件路径
csv_file = '/home/s2331261/Master_Project/6_Feature_Selection/C2_get_active_poses/RDkit_docked_poses_label.csv'
# PDBQT文件夹路径
pdbqt_folder = '/home/s2331261/Master_Project/z1_3p_dataset/dock_results/active_docking_results'
# 输出文件夹路径
output_folder = '/home/s2331261/Master_Project/6_Feature_Selection/C2_get_active_poses/01_split_pdbqt'
# 输出日志文件路径
log_file = '/home/s2331261/Master_Project/6_Feature_Selection/C2_get_active_poses/active_log.txt'


# 创建输出文件夹（如果不存在）
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 读取CSV文件
with open(csv_file, 'r') as file:
    reader = csv.reader(file)
    next(reader)  # 跳过标题行
    for row in reader:
        pose_name = row[0] # 10gs_rdkit_ligand_01
        pdb_code = row[1][-4:] # PDB_10gs; 10gs
        pose_no = int(row[3]) # 1

        pdbqt_file = os.path.join(pdbqt_folder, f'{pdb_code}_rdkit_ligand_out.pdbqt') #1a0q_rdkit_ligand_out.pdbqt
        output_file = os.path.join(output_folder, f'{pdb_code}_pose_{pose_no}.pdbqt')

        # 读取PDBQT文件
        with open(pdbqt_file, 'r') as pdbqt:
            lines = pdbqt.readlines()

        # 查找指定姿势的起始和结束行索引
        start_line = -1
        end_line = -1
        for i, line in enumerate(lines):
            if line.startswith('MODEL') and int(line.split()[1]) == pose_no:
                start_line = i
            elif line.startswith('ENDMDL') and start_line != -1:
                end_line = i
                break

        if start_line != -1 and end_line != -1:
            # 提取选定姿势的行
            selected_lines = lines[start_line:end_line+1]

            # 写入输出文件
            with open(output_file, 'w') as output:
                output.writelines(selected_lines)

            with open(log_file, 'a') as log:
                log.write(f'{pose_name} to {output_file}\n')
        else:
                log.write(f'{pose_name} not found in {pdbqt_file}\n')

print('Extraction complete.')
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))