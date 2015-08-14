#!/bin/bash

for x
do

    ibenchmark.sh $x R1
    ibenchmark.sh $x R2
    ibenchmark.sh $x R3
done
   

exit
