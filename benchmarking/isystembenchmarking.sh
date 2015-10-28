#!/bin/bash

workingdir=./

usage()
{
cat <<EOF
usage:

ISIS: $ISIS

EOF
}

system=`hostname`

#get option
if [ $# -eq 0 ]; then usage; exit 1; fi
while getopts "r:w:" argument; do
    case $argument in
        r ) run=$OPTARG;;
	w ) workflow=$OPTARG;;
        * ) usage
            exit 1;;
    esac
done


#check option
if [[ -z $run ]]; then
    echo "ERROR: Run (-r) required!"
    exit 1
fi

if [[ -z $workflow ]]; then
    echo "ERROR: Run (-w) required!"
    exit 1
fi

processName=${system}${workflow}${run}
validationFile=${workingdir}Aligned/Projects/default/default/sorted.bam.md5

if [ -d "${workingdir}Aligned" ]; then
    echo "ERROR: Aligned diectory all ready exists"
    exit 1
fi

resultsDir=${workingdir}results/
if [ ! -d "$resultDir" ]; then
  mkdir $resultsDir
fi

analysisdir=${workingdir}runs/$processName/
if [ ! -d "$analysisdir" ]; then
    mkdir -p $analysisdir
else
    echo "ERROR: analysis directory already existes"
    exit 1
fi


echo "starting collectl"
collectl -s cmd -f $analysisdir${processName}collectl &
CollectlPid=$!

echo "starting recording cpuMHz"
cpuMHz.sh $analysisdir &
cpuMHzPid=$!

echo "started monerating"


echo "starting iSAAC"
nohup /usr/bin/time -v $workflow.sh > ${analysisdir}${processName}.stdout
echo "finished iSAAC"

echo "killing monerating"
kill $CollectlPid
kill $cpuMHzPid


collectl -p $analysisdir${processName}collectl* -P -f ${analysisdir}collectlplot

gzip -d ${analysisdir}collectlplot*

mv ${analysisdir}collectlplot* ${analysisdir}${processName}.dat

cp $validationFile ${resultsDir}${processName}.val

grep -i "error" ${analysisdir}${processName}.stdout > ${resultsDir}${processName}.error

mv ${workingdir}Aligned ${processName}Aligned

grep "Elapsed (wall clock) time" ${analysisdir}${processName}.stdout | sed "s/\t/$processName: /" >> ${resultsDir}runtimes

collectlTimeDateConverter.py "$processName" "$analysisdir" "$analysisdir"

cp ${analysisdir}${processName}.tsv ${analysisdir}collectl.tsv

gnuisaac.sh ${analysisdir} 

rm ${analysisdir}collectl.tsv

exit
