import pandas as pd
import os
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
from tqdm import tqdm
import timeit



df = pd.read_csv('/home/s2331261/Master_Project/3_Docking/a8_broken_deocy.csv', sep=",")

os.makedirs('a9_files', exist_ok=True)

# 创建一个空的DataFrame来存储错误信息
# error_df = pd.DataFrame(columns=['name', 'smile', 'error'])
# pdb_list_df = pd.DataFrame(columns=['name', 'smile', 'filename'])


def convert_smiles_to_pdb(row):
    Name= row['name']
    Smile= row['smile']
    decoy = Chem.MolFromSmiles(row['smile'])
    name_parts = row['name'].split('_') 
    new_name = f'{name_parts[0]}_decoy_{name_parts[1]}'  
    if decoy is not None:
        decoy_H= Chem.AddHs(decoy)
        ps = AllChem.ETKDGv2()
        ps.useRandomCoords = True
        check = AllChem.EmbedMolecule(decoy_H, ps)
        if check != 0:
            print(f"Embedding failed at index {Name}: {Smile}")
            with open('a9_broken.txt', 'a') as f:
                f.write(f"Embedding failed at index {Name}: {Smile}\n")
            # error_df.loc[len(error_df)] = [Name, Smile, "Embedding failed"]
        else:
            decoy = Chem.RemoveHs(decoy_H)
            Chem.MolToPDBFile(decoy, f'./a9_files/{new_name}.pdb')
            with open('a9_cottected.txt', 'a') as f:
                f.write(f"Corrected decoy at index {Name}: {Smile}\n")
            # pdb_list_df.loc[len(pdb_list_df)] = [Name, Smile,f'{new_name}.pdb']
    else:
        print(f"Invalid SMILES string at index {Name}: {Smile}")
        with open('a9_broken.txt', 'a') as f:
            f.write(f"Invalid SMILES string at index {Name}: {Smile}\n")
        # error_df.loc[len(error_df)] = [Name, Smile, "Invalid SMILES string"]

# use joblib to run the conversions in parallel
start=timeit.default_timer()
Parallel(n_jobs=30)(delayed(convert_smiles_to_pdb)(row) for i, row in tqdm(df.iterrows(), total=df.shape[0]))
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))

# error_df.to_csv('Broken_Decoy_03.csv', index=False)
# pdb_list_df.to_csv('Correct_Decoy_03.csv', index=False)