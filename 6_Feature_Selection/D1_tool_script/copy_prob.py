import shutil
import os

source_folder = "/home/s2331261/Master_Project/6_Feature_Selection/A1_docked_pdbqt"
destination_folder = "/home/s2331261/Master_Project/6_Feature_Selection/A2_prob_pdbqt"

file_names = [
    "2j95_decoy_825",
    "4kwf_decoy_463",
    "1q8u_decoy_308",
    "4zt8_decoy_205",
    "4zt8_decoy_60",
    "4zt8_decoy_342",
    "5nxg_decoy_790",
    "5nxi_decoy_670",
    "5qay_decoy_717",
    "5txy_decoy_187",
    "4tx6_decoy_579",
    "5edb_decoy_518",
    "1gi1_decoy_538",
    "1k4g_decoy_456",
    "1bcu_decoy_239",
    "1bcu_decoy_634",
    "1bgq_decoy_485",
    "1bgq_decoy_42",
    "1bgq_decoy_241",
    "1d6w_decoy_112",
    "1d6w_decoy_220",
    "1dbt_decoy_585",
    "1drv_decoy_583",
    "1e4h_decoy_494",
    "1fcy_decoy_318",
    "1flr_decoy_737",
    "1fvt_decoy_5",
    "2oxd_decoy_379",
    "2oxd_decoy_457",
    "2p95_decoy_848",
    "2za5_decoy_786",
    "2zjf_decoy_12",
    "3d7z_decoy_757",
    "4h3f_decoy_25",
    "4h3g_decoy_128",
    "4h3f_decoy_649",
    "4h3f_decoy_604",
    "4h3f_decoy_649",
    "4h3f_decoy_687",
    "4h3g_decoy_162"
]

for file_name in file_names:
    source_file = os.path.join(source_folder, f"{file_name}_out.pdbqt")
    destination_file = os.path.join(destination_folder, f"{file_name}_out.pdbqt")
    shutil.copyfile(source_file, destination_file)
