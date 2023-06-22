 
import os
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem



# give the smile atoms coordinates and temporary hydrogens for addition of 3D structure
#m1 = Chem.MolFromSmiles("C[S+](=O)(C)C=C(NC1=CC(=C(C=C1)Cl)Cl)[O-]") charge type
#m1 = Chem.MolFromSmiles("CS(C)(C)=O") atmo type

m1 = Chem.MolFromSmiles("CCCCCCCCC(=NCCCCNC(=N)N)N1CCN(C(=O)C(N)C2CCNCC2)CC1")
m1_H= Chem.AddHs(m1)
check = AllChem.EmbedMolecule(m1_H,useRandomCoords=True)
if check != 0:
                print(f"Embedding faild")
#Draw.MolToFile(m1, 'test.png')

# remove the temporary hydrogens from now 3D molecule
m1_H = Chem.RemoveHs(m1_H)
Chem.MolToPDBFile(m1_H, 'm1_H.pdb')

