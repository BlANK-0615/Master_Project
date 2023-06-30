from rdkit import Chem


# generate sdf file,bash command
/home/s2331261/Master_Project/BALLOON1.8.2/balloon -f /home/s2331261/Master_Project/BALLOON1.8.2/MMFF94.mff \
      --nconfs 1 --noGA "CSN1Cc2c3[nH]c(=N)n2C(CC2CCCCC2)=NC1=N3" test2.sdf

# 加载SDF文件
supplier = Chem.SDMolSupplier('test.sdf')

# 读取SDF文件中的第一个（也是唯一的）分子
sdf = supplier[0]
sdf_noH = Chem.RemoveHs(sdf)

# 将分子写入到PDB文件中
pdb = Chem.MolToPDBFile(sdf_noH, 'test01.pdb')

