1bcd_decoy_250.pdbqt
(-5.8465, 0.9810000000000003, 16.8035, 15.705, 16.060000000000002, 15.570999999999998)

 return x_center, y_center, z_center, x_range, y_range, z_range

/home/s2331261/Master_Project/z1_3p_dataset/dock_source/all_receptors_ligands/3d7z/3d7z_receptor.pdbqt

 receptor_file = os.path.join(receptor_folder, pdb_code, f"{pdb_code}_receptor.pdbqt")


<path-to-psovina>/psovina --receptor 1ps3_protein.pdbqt --ligand 1ps3_ligand.pdbqt  --center_x  31.951 --center_y 65.5053 --center_z 7.63888 --size_x  33.452 --size_y 27.612  --size_z  35.136  --num_modes 1 --cpu 8

/home/s2331261/Master_Project/y1_psovina_dock/psovina-2.0/build/linux/release/psovina --receptor 1bcd_receptor.pdbqt --ligand 1bcd_decoy_250.pdbqt  --center_x  -5.8465 --center_y 0.98 --center_z 16.8035 --size_x  15.705 --size_y 16.06  --size_z  15.57 --num_modes 5 


3d7z_decoy_218.pdbqt
(43.578, 30.968000000000004, 33.718, 25.214000000000006, 20.604000000000003, 19.074)/home/s2331261/Master_Project/y1_psovina_dock/psovina-2.0/build/linux/release/psovina --receptor 1bcd_receptor.pdbqt --ligand 1bcd_decoy_250.pdbqt  --center_x  -5.8465 --center_y 0.98 --center_z 16.8035 --size_x  15.705 --size_y 16.06  --size_z  15.57 --num_modes 5 

/home/s2331261/Master_Project/y1_psovina_dock/psovina-2.0/build/linux/release/psovina --receptor 3d7z_receptor.pdbqt --ligand 3d7z_decoy_218.pdbqt --center_x  43.578 --center_y 30.968 --center_z 33.718 --size_x  25.214 --size_y 20.604  --size_z  19.074 --num_modes 5 --out 3d7z_decoy_250_out.pdbqt



% <path-to-AutoDockTools>/prepare_ligand4.py -l 1ps3_ligand.mol2 -o 1ps3_ligand.pdbqt \
    -A 'hydrogens' -U 'nphs_lps_waters'

% <path-to-AutoDockTools>/prepare_receptor4.py -r 1ps3_protein.pdb -o 1ps3_protein.pdbqt \
    -A 'hydrogens' -U 'nphs_lps_waters'


/home/s2331261/Master_Project/3_Docking/SCORCH/utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py \
-r 
                '-l', input_file,
                '-A', 'hydrogens',
                '-o', output_file,
                '-U', 'nphs',


/home/s2331261/Master_Project/BALLOON1.8.2/balloon -f /home/s2331261/Master_Project/BALLOON1.8.2/MMFF94.mff \
    --nconfs 1 --noGA "CC1CCC2C(C)(C)C(O)(C(=O)N=CC=CC(=O)CC3=NC(N[N+](C)=O)=NC3)C3(C)CCC14CC21CCC(O)C1(C)CC43" 5u6g.sdf

obabel test.sdf -O 01.pdb

mol = Chem.MolFromPDBFile('01.pdb')
decoy = Chem.RemoveHs(mol)
