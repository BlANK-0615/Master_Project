# This script is designed to convert PDB files to PDBQT format. 
# It first lists all PDB files from a specified directory.

import os
import subprocess
from joblib import Parallel, delayed
from tqdm import tqdm
import timeit

def convert_pdb_to_pdbqt(pdb_file):
    if pdb_file.endswith('.pdb'):
        # Get the base name of the file (no suffix)
        basename = os.path.splitext(pdb_file)[0]
        # Construct paths for input and output files
        input_file = os.path.join('filepath', pdb_file)
        output_file = os.path.join('filepath', f'{basename}_rdm.pdbqt')
        try:
            # run the conversion command
            command = [
                'filepath/utils/MGLTools-1.5.6/bin/pythonsh',
                'filepath/utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py',
                '-l', input_file,
                '-A', 'hydrogens',
                '-o', output_file,
                '-U', 'nphs',
            ]
            subprocess.run(command, check=True)
        except:
            print(f'ERROR: {pdb_file}')
            with open('/filepath/qdbqt_error.txt', 'a') as f:
                f.write(f"Error pdb file at index {pdb_file}\n")


# list pdb files
pdb_files = os.listdir('path_to_pdb_files')

# Convert all pdb files in parallel
# The n_jobs parameter can be adjusted according to your actual situation, such as your number of physical cores or logical cores
start=timeit.default_timer()
Parallel(n_jobs=-1)(delayed(convert_pdb_to_pdbqt)(pdb_file) for pdb_file in tqdm(pdb_files, total=len(pdb_files)))
end=timeit.default_timer()
print('Running time: %s Seconds'%(end-start))

