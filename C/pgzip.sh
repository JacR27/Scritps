#!/bin/bash
tiles=24
swaths=2
surface=2
cycles=302
lane=5
extention=bcl.gz.filtered


for i in `seq -f "%02g" 1 $tiles`; 
do
    for k in `seq 1 $swaths`; 
    do
	for j in `seq 2 $surface`;
	do
	    for n in `seq 1 $cycles`;
	    do
		gzip "C${n}.1/s_${lane}_${j}${k}${i}.$extention" &
		pid=$!
	    done
	    wait $pid
	done
    done
done
exit
