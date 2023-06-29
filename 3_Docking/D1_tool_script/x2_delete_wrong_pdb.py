import os

folder_path = '/home/s2331261/Master_Project/3_Docking/pdb_files'  # 文件夹路径

# 获取文件夹内所有文件名
file_names = os.listdir(folder_path)

# 遍历文件名，删除问题文件
for file_name in file_names:
    if '_decoy_decoy.pdb' in file_name:
        file_path = os.path.join(folder_path, file_name)
        os.remove(file_path)
        print(f"Deleted file: {file_name}")

# Deleted file: 3u90_decoy_decoy.pdb
# Deleted file: 6gon_decoy_decoy.pdb
# Deleted file: 3w9k_decoy_decoy.pdb
# Deleted file: 4qf8_decoy_decoy.pdb
# Deleted file: 4pox_decoy_decoy.pdb
# Deleted file: 4d3h_decoy_decoy.pdb
# Deleted file: 6h8s_decoy_decoy.pdb
# Deleted file: 5eqy_decoy_decoy.pdb
# Deleted file: 5org_decoy_decoy.pdb
# Deleted file: 4q3u_decoy_decoy.pdb
# Deleted file: 3v2q_decoy_decoy.pdb
# Deleted file: 5ot8_decoy_decoy.pdb
# Deleted file: 4pp0_decoy_decoy.pdb
# Deleted file: 1hmt_decoy_decoy.pdb
# Deleted file: 2gv6_decoy_decoy.pdb
# Deleted file: 3uex_decoy_decoy.pdb
# Deleted file: 2aco_decoy_decoy.pdb
# Deleted file: 1g74_decoy_decoy.pdb
# Deleted file: 3usx_decoy_decoy.pdb
# Deleted file: 4nja_decoy_decoy.pdb
# Deleted file: 4nj9_decoy_decoy.pdb
# Deleted file: 4l40_decoy_decoy.pdb
# Deleted file: 4tkj_decoy_decoy.pdb
# Deleted file: 3kwa_decoy_decoy.pdb
# Deleted file: 4rww_decoy_decoy.pdb
# Deleted file: 3wvm_decoy_decoy.pdb
# Deleted file: 5ito_decoy_decoy.pdb
# Deleted file: 4pow_decoy_decoy.pdb
# Deleted file: 4tkb_decoy_decoy.pdb
# Deleted file: 3qlm_decoy_decoy.pdb
# Deleted file: 4gny_decoy_decoy.pdb
# Deleted file: 1pot_decoy_decoy.pdb
# Deleted file: 1vyf_decoy_decoy.pdb
# Deleted file: 6ge7_decoy_decoy.pdb
# Deleted file: 1hms_decoy_decoy.pdb
# Deleted file: 2hnx_decoy_decoy.pdb
# Deleted file: 5itp_decoy_decoy.pdb
# Deleted file: 6gf9_decoy_decoy.pdb
# Deleted file: 4fnn_decoy_decoy.pdb
# Deleted file: 1w7g_decoy_decoy.pdb
# Deleted file: 2ifb_decoy_decoy.pdb
# Deleted file: 1gt4_decoy_decoy.pdb
# Deleted file: 3uev_decoy_decoy.pdb
# Deleted file: 4rse_decoy_decoy.pdb
# Deleted file: 3v2p_decoy_decoy.pdb
# Deleted file: 1hmr_decoy_decoy.pdb
# Deleted file: 4wk1_decoy_decoy.pdb
# Deleted file: 3uew_decoy_decoy.pdb
# Deleted file: 5uxf_decoy_decoy.pdb
# Deleted file: 5f1v_decoy_decoy.pdb
# Deleted file: 5ijr_decoy_decoy.pdb
# Deleted file: 6gfs_decoy_decoy.pdb