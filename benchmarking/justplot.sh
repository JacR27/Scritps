#!/bin/bash
PATH=$PATH:/mnt/nfs/Scritps/benchmarking
export PATH

analysisdir=./


if [ $# -eq 0 ]; then usage; exit 1; fi
while getopts "d:" argument; do
    case $argument in
        d ) processName=$OPTARG;;
        * ) usage
            exit 1;;
    esac
done

echo $processName
echo $analysisdir

collectl -p $analysisdir${processName}collectl* -P -f ${analysisdir}collectlplot

echo "converting to plotable file"

gzip -d ${analysisdir}collectlplot*

echo "unziping plotable file"

mv ${analysisdir}collectlplot* ${analysisdir}${processName}.dat

echo "renameing plotable file"

collectlTimeDateConverter.py "$processName" "$analysisdir" "$analysisdir"

echo "converting dates"

grep "Elapsed time" ${analysisdir}WorkflowLog.txt | extractIsisSteps.py > ${analysisdir}${processName}.stp

echo "getting subprocessing"

cp ${analysisdir}${processName}.tsv ${analysisdir}collectl.tsv

echo "renameing file"

gnuisaac.sh ${analysisdir}

echo "plotting"

rm ${analysisdir}collectl.tsv

echo "removeing temp file"

exit
