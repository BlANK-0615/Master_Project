import sys
sys.path.append("../")
import binana
ligand, receptor = binana.load_ligand_receptor.from_files("/home/s2331261/Master_Project/6_Feature_Selection/z1_binana_test/6ekq/ligands/6ekq_decoy_1_ligand.pdbqt", "/home/s2331261/Master_Project/6_Feature_Selection/z1_binana_test/6ekq/6ekq_receptor.pdbqt")

ligand, receptor
# This script could be used under /home/s2331261/Master_Project/6_Feature_Selection/binana/python