#!/bin/bash

system=`hostname`
workflow=$1
run=$2
processName=${system}${workflow}${run}
workingdir=./
analysisdir=$workingdir/runs/$processName/
command=/illumina/development/iSAAC/iSAAC-03.15.07.15/bin/isaac-align


echo "this script will run iSAAC benchmarking and monitor system using collectl"

mkdir $analysisdir

collectl -s cmd -f $workingdir$analysisdir$processName &
CollectlPid=$!

echo collectl started

nohup /usr/bin/time -v $command -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format bcl-gz -m 59  --tiles s_4 -o ${workingdir}${analysisdir}/Align > ${workingdir}${analysisdir}${processName}.stdout

echo killing collectl

kill $CollectlPid

collectl -p $workingdir$analysisdir${processName}* -P -f ${workingdir}${analysisdir}plot*

gzip -d ${workingdir}${analysisdir}plot*

mv ${workingdir}${analysisdir}plot* ${workingdir}${analysisdir}${processName}.dat
