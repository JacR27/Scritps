#!/bin/bash

for x
do

    /home/sbsuser/Scritps/benchmarking/isystembenchmarking.sh $x R1
    /home/sbsuser/Scritps/benchmarking/isystembenchmarking.sh $x R2
    /home/sbsuser/Scritps/benchmarking/isystembenchmarking.sh $x R3

done
   

exit
