#!/bin/bash

mkdir -p pdbqt_files

for file in pdb_files/*.pdb; do
    base=$(basename "$file" .pdb)
    base_parts=(${base//_/ })  # Split the base name at underscore
    new_base="${base_parts[0]}_decoy_${base_parts[1]}"  # Combine with _decoy_
    ./utils/MGLTools-1.5.6/bin/pythonsh ./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l "$file" -A hydrogens -o "pdbqt_files/${new_base}.pdbqt" -U nphs
done


#!/bin/bash

mkdir -p pdbqt_files

files=(pdb_files/*.pdb)
total=${#files[@]}

for i in "${!files[@]}"; do
    file=${files[$i]}
    base=$(basename "$file" .pdb)
    base_parts=(${base//_/ })  # Split the base name at underscore
    new_base="${base_parts[0]}_decoy_${base_parts[1]}"  # Combine with _decoy_
    ./utils/MGLTools-1.5.6/bin/pythonsh ./utils/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l "$file" -A hydrogens -o "pdbqt_files/${new_base}.pdbqt" -U nphs

    echo "Processed: $(($i+1))/$total"
done
