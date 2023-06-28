#!/usr/people/douglas/programs/ksh.exe
export LD_LIBRARY_PATH=/usr/people/douglas/programs/boost_1_59_0/lib


function vina_prep {
if [[ ! -e vina.dpf ]]; then
  error="Vina.dpf not found; please run autoAD.ksh -j to generate it."
  errhan $0
fi
centers=($(grep gridcenter vina.dpf))
npts=($(grep npts vina.dpf))
receptor=$(grep .pdbqt vina.dpf)
max_params=$(grep exhaustiveness vina.dpf)
max_params="--exhaustiveness=32"
if (( pose == 1 )); then
  prepose="--local_only"
fi
spacing=0.2
rm -f ${arglist[0]%.pdbqt}_out.dlg
}

function flags {
print -- "--receptor ${dir}_protein_H.pdbqt --ligand ${dir}_ligand.pdbqt --center_x ${centers[1]} --center_y ${centers[2]} --center_z ${centers[3]} --size_x $(((${npts[1]}+1)*$spacing)) --size_y $(((${npts[2]}+1)*$spacing)) --size_z $(((${npts[3]}+1)*$spacing)) "$max_params" $prepose "
}


rm -f RMSD_results_*max.txt
dirs=$(ls)
print dirs are $dirs
for dir in $dirs; do
  if [[ -d $dir ]]; then
    cd $dir
#####################################################################WATCH THIS LINE
#    print "Deleting all dlg files"
#    rm -rf *.dlg
    print "dir is $dir"


#    print $dir $(~/scripts/rmsd_min.ksh ${dir}_ligand.pdbqt ${dir}_ligand_out_vina112.pdbqt) >> ../RMSD_results_vina112.txt

    print $dir $(~/scripts/rmsd_min.ksh ${dir}_ligand.pdbqt ${dir}_ligand_out_psovina2max.pdbqt) >> ../RMSD_results_psovina2max.txt
    print $dir $(~/scripts/rmsd_min.ksh ${dir}_ligand.pdbqt ${dir}_ligand_out_gwovinamax.pdbqt) >> ../RMSD_results_gwovinamax.txt

#    print $dir $(~/scripts/rmsd_min.ksh ${dir}_ligand.pdbqt ${dir}_ligand_out_quickvina2.pdbqt) >> ../RMSD_results_quickvina2.txt
#    print $dir $(~/scripts/rmsd_min.ksh ${dir}_ligand.pdbqt ${dir}_ligand_out_quickvinaw.pdbqt) >> ../RMSD_results_quickvinaw.txt


#    ~/scripts/docking/autoAD_test.ksh -qpj ${dir}_ligand.mol2 ${dir}_protein.pdb ${dir}_ligand.mol2 
#    vina_prep
#    print -- "$(flags)"
#    time -o time_gwovinamax.txt /usr/people/douglas/programs/gwovina-1.0/build/linux/release/gwovina $(print -- "$(flags)") | tee -a ${dir}_out_gwovina.dlg; mv ${dir}_ligand_out.pdbqt ${dir}_ligand_out_gwovinamax.pdbqt
#    time -o time_vina112.txt  /usr/people/douglas/programs/autodock_vina_1_1_2_linux_x86/bin/vina $(print -- "$(flags)") | tee -a ${dir}_out_vina112.dlg; mv ${dir}_ligand_out.pdbqt ${dir}_ligand_out_vina112.pdbqt 
#    time -o time_psovina2max.txt  /usr/people/douglas/programs/psovina2 $(print -- "$(flags)") | tee -a ${dir}_out_psovina2.dlg; mv ${dir}_ligand_out.pdbqt ${dir}_ligand_out_psovina2max.pdbqt
#    time -o time_qvina2.1.txt  /usr/people/douglas/programs/QuickVina/qvina2.1 $(print -- "$(flags)") | tee -a ${dir}_out_quickvina2.dlg; mv ${dir}_ligand_out.pdbqt ${dir}_ligand_out_quickvina2.pdbqt
#    time -o time_qvinaw.txt  /usr/people/douglas/programs/QuickVina/qvina-w $(print -- "$(flags)") | tee -a ${dir}_out_quickvinaw.dlg; mv ${dir}_ligand_out.pdbqt ${dir}_ligand_out_quickvinaw.pdbqt

    cd ..
  fi
done
