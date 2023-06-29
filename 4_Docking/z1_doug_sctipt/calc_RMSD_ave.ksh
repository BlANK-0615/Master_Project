#!/usr/people/douglas/programs/ksh.exe


for file in $( find CASF-2016_QuickVina2_W_PSOVina2_GWOVina_Vina112/coreset/ -name "RMSD_results*max_cleaned.txt"); do
total=0
c=0
success=0
for val in $(awk '{print $2}' $file); do
  ((total+=val))
  ((c++))
  if ((val <= 2.00 )); then
    ((success++))
  fi
done


average=$((total/${c}.0))
successrate=$(( (success/${c}.0)*100))
printf "$file average is %.2f; success rate is %.2f\%\n" $average $successrate
done
