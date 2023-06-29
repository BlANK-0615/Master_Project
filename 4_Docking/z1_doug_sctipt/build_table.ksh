#!/usr/people/douglas/programs/ksh.exe



function my_times {
for dir in $(ls -1d -- */ | grep -v Vina); do
  total=0
  seconds=0
  minutes=0
  hours=0
  second=0
  minute=0
  hour=0
#  print dir is $dir
  seconds=$(grep $dir alltimes.txt | awk '{print $3}' | sed "s/elapsed//" | sed "s/\..*//" | sed "s/.*\://"  )
  minutes=$(grep $dir alltimes.txt | awk '{print $3}' | sed "s/elapsed//" | sed "s/\..*//" | sed "s/.*\:\(.*\)\:.*/\1/" | sed "s/\:.*//"  )
  hours=$(grep $dir alltimes.txt | awk '{print $3}' | sed "s/elapsed//" | sed "s/\..*//" | sed 's/:[[:digit:]]*$//' | sed 's/\(.*[[:digit:]]\:\)[[:digit:]]*/\1/' | sed 's/[[:digit:]]*[^\:]$/0/' | sed "s/\://")

#  sleep 1
#  print for $dir total is $total 
  for second in $seconds; do
#    print second is $second
    ((total+=second))
  done
  for minute in $minutes; do
#    print minute is $minute
#    print "((total+=(${minute-0}*60)))"
    ((total+=(${minute-0}*60)))
  done
  for hour in $hours; do
#    print hour is $hour
    ((total+=(${hour-0}*3600)))
  done
  print for $dir total is $total

done

}

function my_times_vina {
for dir in $(ls -1d -- */ | grep Vina); do
  for file in vina112 qvinaw qvina2 psovina2 gwovina; do
  total=0
  seconds=0
  minutes=0
  hours=0
  second=0
  minute=0
  hour=0
#  print dir is $dir
  seconds=$(grep $dir  alltimes.txt | grep $file | awk '{print $3}' | sed "s/elapsed//" | sed "s/\..*//" | sed "s/.*\://"  )
  minutes=$(grep $dir alltimes.txt | grep $file | awk '{print $3}' | sed "s/elapsed//" | sed "s/\..*//" | sed "s/.*\:\(.*\)\:.*/\
\1/" | sed "s/\:.*//"  )
  hours=$(grep $dir alltimes.txt | grep $file |  awk '{print $3}' | sed "s/elapsed//" | sed "s/\..*//" | sed 's/:[[:digit:]]*$//'\
 | sed 's/\(.*[[:digit:]]\:\)[[:digit:]]*/\1/' | sed 's/[[:digit:]]*[^\:]$/0/' | sed "s/\://")

#  sleep 1
#  print for $dir $file total is $total 
  for second in $seconds; do
#    print second is $second
    ((total+=second))
  done
  for minute in $minutes; do
#    print minute is $minute
#    print "((total+=(${minute-0}*60)))"
    ((total+=(${minute-0}*60)))
  done
  for hour in $hours; do
#    print hour is $hour
    ((total+=(${hour-0}*3600)))
  done
  print for $dir $file total is $total

done
done
}
my_times
my_times_vina
