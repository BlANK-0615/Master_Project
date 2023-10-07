# This script reads a CSV file containing molecular names 
# and their respective SMILES strings. For each molecule, 
# it attempts to convert the SMILES to a 3D structure and then save as a PDB file.


import pandas as pd
import os
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
from tqdm import tqdm
import timeit

# load the CSV file containing the SMILES strings
df = pd.read_csv('path/to/file', sep=",")
os.makedirs('', exist_ok=True)

# create a dataframe to store the results
def convert_smiles_to_pdb(row):
    Name= row['name']
    Smile= row['smile']
    decoy = Chem.MolFromSmiles(row['smile'])
    name_parts = row['name'].split('_') 
    new_name = f'{name_parts[0]}_decoy_{name_parts[1]}'  
    if decoy is not None:
        # embed the molecule in 3D space
        decoy_H= Chem.AddHs(decoy)
        ps = AllChem.ETKDGv2()
        ps.useRandomCoords = True
        check = AllChem.EmbedMolecule(decoy_H, ps)
        # if embedding fails, write the name and SMILES string to a text file
        if check != 0:
            print(f"Embedding failed at index {Name}: {Smile}")
            with open('broken.txt', 'a') as f:
                f.write(f"Embedding failed at index {Name}: {Smile}\n")
            # error_df.loc[len(error_df)] = [Name, Smile, "Embedding failed"]
        else: 
            # write the decoy to a PDB file
            decoy = Chem.RemoveHs(decoy_H)
            Chem.MolToPDBFile(decoy, f'./filepath/{new_name}.pdb')
            with open('correct.txt', 'a') as f:
                f.write(f"Corrected decoy at index {Name}: {Smile}\n")
            # pdb_list_df.loc[len(pdb_list_df)] = [Name, Smile,f'{new_name}.pdb']
    else:
        print(f"Invalid SMILES string at index {Name}: {Smile}")
        with open('broken.txt', 'a') as f:
            f.write(f"Invalid SMILES string at index {Name}: {Smile}\n")
        # error_df.loc[len(error_df)] = [Name, Smile, "Invalid SMILES string"]

# use joblib to run the conversions in parallel
start=timeit.default_timer()
Parallel(n_jobs=30)(delayed(convert_smiles_to_pdb)(row) for i, row in tqdm(df.iterrows(), total=df.shape[0]))
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))

