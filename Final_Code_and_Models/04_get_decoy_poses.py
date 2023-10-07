# The main purpose of this script is to randomly select a pose
# for each PDBQT file from multiple poses and record the selected pose to a CSV file.

import random
import pandas as pd
import os
import shutil
import timeit

# get the docked poses from the output pdbqt files 
def parse_pdbqt(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    poses = []
    pose = []
    for line in lines:
        if line.startswith("MODEL"):
            pose = [line.strip()]  # Add the MODEL line to the pose
        elif line.startswith("ENDMDL"):
            pose.append(line.strip())  # Add the ENDMDL line to the pose
            poses.append(pose)
            pose = []
        else:
            pose.append(line.strip())
            
    return poses

# randomly select a pose from the list of poses and write it to a file
def random_select_and_record(poses, file_path, csv_path, output_folder):
    results = {}
  
    if len(poses) == 0:
       # Handle invalid files and record error messages
        error_message = "No valid poses found in the file"
        results[file_path] = error_message
        with open("Error_poses.txt", 'a') as error_file:
            error_file.write(f"{file_path}: {error_message}\n")
    
    # randomly select a pose
    selected_pose_index = random.randint(0, len(poses) - 1)
    selected_pose = poses[selected_pose_index]
    
    # get the pose name and number
    pose_name = os.path.basename(file_path).replace("_out.pdbqt", "")
    pose_no = selected_pose_index + 1
    
    # write the selected pose to a file
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


# Process all pdbqt files in the specified folder, 
# do a random selection of poses for each file, 
# and log the results to a CSV file.
# Also output the selected poses as separate pdbqt files to another directory.
def process_files(directory, output_folder, output_csv):
    # make sure the output folder exists
    os.makedirs(output_folder, exist_ok=True)
    # obtain all file names in the specified directory
    file_names = os.listdir(directory)
    # filter out all pdbqt files
    pdbqt_files = [f for f in file_names if f.endswith('.pdbqt')]

    for file_name in pdbqt_files:
        # obtain the full path of the file
        file_path = os.path.join(directory, file_name)
        # parse the file and obtain all poses
        poses = parse_pdbqt(file_path)
        # do a random selection of poses and record the results
        random_select_and_record(poses, file_path, output_csv, output_folder)

# run the process_files function
start=timeit.default_timer()
process_files("path_to_pdbqt_files",
                "path_to_output_folder", 
                "path_to_output_csv_file")

end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))
