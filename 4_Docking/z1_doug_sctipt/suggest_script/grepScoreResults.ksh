#!/usr/people/douglas/programs/ksh.exe

rm -f *_scores_max.txt
set -A score
for dir in $(< dirlist.txt); do
  print "dir is $dir"
  for file in gwovinamax psovina2max; do
    score[$file]=$(head -2 ${dir}/${dir}_ligand_out_${file}.pdbqt | tail -1 | awk '{print $4}')
    print -- $dir ${score[$file]} >> ${file}_scores_max.txt
  done
done

for file in gwovinamax psovina2max; do
   sort -n --key=2,2 ${file}_scores_max.txt > ${file}_scores_max_sorted.txt
done 
