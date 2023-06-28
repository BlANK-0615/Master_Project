#!/usr/people/douglas/programs/ksh.exe


set -A score
rm -f *_scores.txt
for dir in $(< dirlist.txt); do
  print "dir is $dir"
  for file in CASF-2016_AutoDockAMBERcharges CASF-2016_AutoDockAMBERcharges_maxparams CASF-2016_AutoDockPARSEcharges CASF-2016_AutoDockPARSEcharges_maxparams; do
    score[${file//\-}]=$(grep "Estimated Free Energy of Binding" ${file}/coreset/${dir}/${dir}_ligand1_largestC.pdbqt | awk '{print $8}')
    print -- $dir ${score[${file//\-}]} >> ${file//\-}_scores.txt
  done
done


for file in CASF-2016_AutoDockAMBERcharges CASF-2016_AutoDockAMBERcharges_maxparams CASF-2016_AutoDockPARSEcharge\
s CASF-2016_AutoDockPARSEcharges_maxparams; do
  sort -n --key=2,2 ${file//\-}_scores.txt | grep -v 3ag9 | grep -v 3utu | grep -v  3zso | grep -v  4ciw | grep -v  4llx | grep -v  4mme | grep -v 4twp > ${file//\-}_scores_sorted2.txt
  printf "${file//\-}: "
  ~/scripts/docking/spearmans.ksh  ${file//\-}_scores_sorted2.txt realaffinities_sorted_nofails.txt | tail -1
done
