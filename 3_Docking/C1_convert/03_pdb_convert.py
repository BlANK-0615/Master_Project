import pandas as pd
import os
from tqdm import tqdm
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem

df = pd.read_csv('/home/s2331261/Master_Project/3_Docking/A5_4ojz/4ojz.csv', sep=",")

os.makedirs('/home/s2331261/Master_Project/3_Docking/A5_4ojz/4ojz_pdb', exist_ok=True)

with tqdm(total=len(df['smile'])) as pbar:
    for index, row in df.iterrows():
        # Original molecule name: 3ljg_0
        Name= row['name']  # Result: 3ljg_0
        Smile= row['smile']  # Result: Clc1ccc2c(c1)nccc2Sc1nnc(s1)NC(=O)c1cccs1
        # name_parts = row['name'].split('_')  # Result: ['3ljg', '0']
        # new_name = f'{name_parts[0]}_decoy_{name_parts[1]}'  # Result: 3ljg_decoy_0
        # new_name = f'{Name}_deep' 
        decoy = Chem.MolFromSmiles(row['smile'])
        if decoy is not None:
            # Add temporary hydrogens & remove before saving to pdb file
            decoy_H= Chem.AddHs(decoy)
            check = AllChem.EmbedMolecule(decoy_H, randomSeed=0xf00d)
            if check != 0:
                print(f"Embedding failed at index {Name}: {Smile}")
                with open('/home/s2331261/Master_Project/3_Docking/A5_4ojz/Broken_Decoy.txt', 'a') as f:
                    f.write(f"Embedding failed at index {Name}: {Smile}\n")
            else:
                decoy = Chem.RemoveHs(decoy_H)
                Chem.MolToPDBFile(decoy, f'/home/s2331261/Master_Project/3_Docking/A5_4ojz/4ojz_pdb/{Name}.pdb')
        else:
            print(f"Invalid SMILES string at index {Name}: {Smile}")
            with open('/home/s2331261/Master_Project/3_Docking/A5_4ojz/Broken_Decoy.txt', 'a') as f:
                f.write(f"Invalid SMILES string at index {Name}: {Smile}\n")
        pbar.update(1)

