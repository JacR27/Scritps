#!/bin/bash
tiles=24
swaths=1
surface=1
cycles=302
lane=5
extention=bcl.gz.rm


for i in `seq -f "%02g" 1 $tiles`; 
do
    for k in `seq 1 $swaths`; 
    do
	for j in `seq 1 $surface`;
	do
	    for n in `seq 1 $cycles`;
	    do
		gzip "C${n}.1/s_${lane}_${j}${k}${i}.$extention" &
	    done
	    wait $pid
	done
    done
done
exit
