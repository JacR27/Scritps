#!/bin/bash 

analysisdir=$1

while [ 1 ]; do
    sleep 1
    dt=`date`
    totalMHz=`grep 'MHz' /proc/cpuinfo | python cpuMHz.py`
    echo "dt: $totalMHz" | sed "s/dt/$dt/" >>${analysisdir}MHz
    
    date >> ${analysisdir}MHzBD
    grep 'MHz' /proc/cpuinfo >> ${analysisdir}MHzBD
done

echo done
