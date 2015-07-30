#!/bin/bash

system=`hostname`
workflow=$1
run=$2
processName=${system}${workflow}${run}
workingdir=./
analysisdir=$workingdir/runs/$processName/
command=/illumina/development/Isis/2.5.55.16.NorthStar/Isis
sampleSheetdir=~/sripts/benchmarking/SampleSheets/

echo "this script will run iSAAC benchmarking and monitor system using collectl"

cp ${sampleSheetdir}SampleSheet.csv.$workflow $workingdir/SampleSheet.csv

mkdir $analysisdir

collectl -s cmd -f $workingdir$analysisdir$processName &
CollectlPid=$!

echo collectl started

nohup /usr/bin/time -v $command -r $workingdir -a $analysisdir> ${workingdir}${analysisdir}${processName}.stdout

echo killing collectl

kill $CollectlPid

collectl -p $workingdir$analysisdir${processName}* -P -f ${workingdir}${analysisdir}plot*

gzip -d ${workingdir}${analysisdir}plot*

mv ${workingdir}${analysisdir}plot* ${workingdir}${analysisdir}${processName}.dat
