#!/bin/bash

system=`hostname`
workflow=$1
run=$2
processName=${system}${workflow}${run}
workingdir=./
analysisdir=${workingdir}runs/$processName/
validationFile=${workingdir}Aligned/Projects/default/default/sorted.bam.md5
resultsDir=${workingdir}results/

echo "this script will run iSAAC benchmarking and monitor system using collectl"

mkdir $analysisdir

collectl -s CMDZ -i5:50 -f $analysisdir$processName &
CollectlPid=$!

echo collectl started

nohup /usr/bin/time -v $workflow.sh > ${analysisdir}${processName}.stdout

echo killing collectl

kill $CollectlPid

mv ${workingdir}Aligned ${processName}Aligned

rm ${workingdir}Temp/{bin-*,gnuplot-*}


exit
