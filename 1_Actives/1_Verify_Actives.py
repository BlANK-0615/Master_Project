import pandas as pd

df = pd.read_csv('/home/s2331261/Master_Project/1_Actives/Active_Ligand_Smiles.csv', dtype={'PDBCode': str})

for index, row in df.iterrows():
    pdb_code = row['PDBCode']

    if not pd.Series(pdb_code).str.match('\d\w{3}').all():
        print(f'Row {index} has an invalid PDB code: {pdb_code}')

