output = subprocess.check_output(f'python {binana_path} -receptor {receptor_path} -ligand {ligand_path} > {destination_path}{ligand}.txt


python3 /home/s2331261/Master_Project/z2_3p_software/Other/binana.py \
    -receptor /home/s2331261/Master_Project/6_Feature_Selection/6ekq/6ekq_receptor.pdbqt \
    -ligand /home/s2331261/Master_Project/6_Feature_Selection/6ekq/ligands/6ekq_decoy_1_ligand.pdbqt \
    > 6ekq_B_test.txt

# run BINANA2.1
python3 /home/s2331261/Master_Project/6_Feature_Selection/binana/python/run_binana.py \
    -receptor /home/s2331261/Master_Project/6_Feature_Selection/6ekq/6ekq_receptor.pdbqt \
    -ligand /home/s2331261/Master_Project/6_Feature_Selection/6ekq_decoy_193_out.pdbqt \
    -output_dir /home/s2331261/Master_Project/6_Feature_Selection/6ekq_Toco_test > 6ekq_Toco_errors.txt


python3 /home/milesm/Dissertation/Code/featureExtraction/calculate_ecifs.py 
-pdbqts /home/milesm/Dissertation/Data/PDBQT/Non_redundant/Druglike_filtered_binding_data_with_decoys/ 
-pdbs /home/milesm/Dissertation/Data/PDBQT/Non_redundant/pdb_copies_of_druglike/ 
-out /home/milesm/Dissertation/Data/csv/ecif_data.csv

python3 /home/s2331261/Master_Project/6_Feature_Selection/03_calculate_ecifs.py \
-pdbqts /home/s2331261/Master_Project/6_Feature_Selection/6ekq \
-pdbs /home/s2331261/Master_Project/6_Feature_Selection/6ekq \
-out /home/s2331261/Master_Project/6_Feature_Selection/6ekq_test.csv

scp -r /home/s2331261/Master_Project/3_Docking/A6_B_Atom s2331261@executor.bch.ed.ac.uk:/home/s2331261/Msater_Project/dock_job/