#!/usr/people/douglas/programs/ksh.exe


function RMSD {
filelist=$(find . -name "RMSD_results*")
for dir in $(< dirlist.txt); do
  for file in $filelist; do
    printf "dir is $dir, file is $file;  "
    grep $dir $file
  done

done
}
