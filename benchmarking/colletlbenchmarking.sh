#!/bin/bash

usage()
{
cat <<EOF
usage:

ISIS: $ISIS

EOF
}


system="S4"
workingdir=./
continueFlag=1
#get option
if [ $# -eq 0 ]; then usage; exit 1; fi
while getopts "r:w:c" argument; do
    case $argument in
        r ) run=$OPTARG;;
        w ) workflow=$OPTARG;;
        c ) continueFlag=0;;
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

resultsDir=${workingdir}results/
if [ ! -d "$resultsDir" ]; then
  mkdir $resultsDir
fi

analysisdir=${workingdir}$processName/
if [ $continueFlag = "1" ]; then
    if [ ! -d "$analysisdir" ]; then
        mkdir -p $analysisdir
    else
        echo "ERROR: analysis directory already existes"
        exit 1
    fi
fi




collectl -s cmd -f $analysisdir${processName}collectl &
CollectlPid=$!
sleep 30
kill $CollectlPid

collectl -p $analysisdir${processName}collectl* -P -f ${analysisdir}collectlplot

gzip -d ${analysisdir}collectlplot*

mv ${analysisdir}collectlplot* ${analysisdir}${processName}.dat

exit

