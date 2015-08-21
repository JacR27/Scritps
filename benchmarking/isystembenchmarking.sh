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

collectl -s cmd -f $analysisdir$processName &
CollectlPid=$!

echo collectl started

nohup /usr/bin/time -v /home/sbsuser/Scritps/benchmarking/$workflow.sh > ${analysisdir}${processName}.stdout

echo killing collectl

kill $CollectlPid

collectl -p $analysisdir${processName}* -P -f ${analysisdir}plot*

gzip -d ${analysisdir}plot*

mv ${analysisdir}plot* ${analysisdir}${processName}.dat

cp $validationFile ${resultsDir}${processName}.val

grep -i "error" ${analysisdir}${processName}.stdout > ${resultsDir}${processName}.error

mv ${workingdir}Aligned ${processName}Aligned

grep "Elapsed (wall clock) time" ${analysisdir}${processName}.stdout | sed "s/\t/$processName: /" >> ${resultsDir}runtimes

rm ${workingdir}Temp/{bin-*,gnuplot-*}

scriptsDir=/home/sbsuser/Scritps/benchmarking/

grep "Elapsed time" ${analysisdir}${processName}.stdout | python ${scriptsDir}extractIsisSteps.py > ${analysisdir}$processName

python ${scriptsDir}collectlTimeDateConverter.py "$processName" "$analysisdir" "$analysisdir"

cp ${analysisdir}${processName}.tsv ${analysisdir}collectl.tsv

/home/sbsuser/Scripts/benchmarking/gnuisaac.sh ${analysisdir} 

rm ${analysisdir}collectl.tsv

exit
