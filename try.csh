#!/bin/csh

#set arr = set `/usr/bin/seq 1 10`
set arr = (0 1 2 3 4 5 6 7 8 9)
set i = 1
set ei = 11
while ($i < $ei)
#echo ${arr[$i]}
set yy = `/usr/bin/printf "%04d" $i`
set yy = `/usr/bin/expr $yy - 1`
set arr[$i]=$yy
echo ${arr[$i]}
set i = `/usr/bin/expr $i + 1`


end


