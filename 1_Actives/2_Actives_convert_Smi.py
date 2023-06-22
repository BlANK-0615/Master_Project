import numpy as np
import pandas as pd


# read the SMILES file
df = pd.read_csv('Active_Ligand_Smiles.csv',dtype={'PDBCode':str})

# check the broken PDBcode
df[df['DatabaseSmile']=='O[C@@H]([C@H]([C@@H]([C@@H](CC(=O)C(=O)O)O)O)O)COP(O)(O)O']
#df[df['PDBCode']=='OC(=O)CCC(=O)C(=O)O']

# get the neeed column
df_RDK=df.drop(labels=['RDKitCanonSmile'],axis=1) 
df_DATA=df.drop(labels=['DatabaseSmile'],axis=1) 

# change the column order, make the format for TocoDecoy generation 
df_RDK[['PDBCode','DatabaseSmile']]=df_RDK[['DatabaseSmile','PDBCode']]
df_DATA[['PDBCode','RDKitCanonSmile']]=df_DATA[['RDKitCanonSmile','PDBCode']]

# check the format
df_RDK
df_DATA

# check the broken PDBcode
df_RDK[df_RDK['PDBCode']=='O[C@@H]([C@H]([C@@H]([C@@H](CC(=O)C(=O)O)O)O)O)COP(O)(O)O']
df_DATA[df_DATA['RDKitCanonSmile']=='3e12']

# save the file to smi format to geanerate the TocoDecoy
df_RDK.to_csv('RDK_Smiles.smi',sep=' ',index=0,header=0)
df_DATA.to_csv('Data_Smiles.smi',sep=' ',index=0,header=0)
