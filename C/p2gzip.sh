#!/bin/bash
cycles=302

for n in `seq 151 $cycles`;
do
    gunzip C${n}.1/* &
    pid=$!
done
wait $pid


exit
