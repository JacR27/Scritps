#!/bin/bash
PATH=$PATH:/mnt/nfs/Scritps/benchmarking
export PATH
workingdir=/mnt/raid/150509_HSX213_0053_AH5WWCCCXX/
usage()
{
cat <<EOF
usage:

ISIS: $ISIS

EOF
}



system="S4"
con=1
#get option
if [ $# -eq 0 ]; then usage; exit 1; fi
while getopts "r:w:p:c" argument; do
    case $argument in
        r ) run=$OPTARG;;
	w ) workflow=$OPTARG;;
	p ) processor=$OPTARG;;
	c ) con=0;;
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

if [[ -z $processor ]]; then
    echo "ERROR: Run (-p) required!"
    exit 1
fi

processName=${system}${workflow}${run}


resultsDir=${workingdir}results/
if [ ! -d "$resultsDir" ]; then
  mkdir $resultsDir
fi

analysisdir=${workingdir}$processName/
echo $con
if [ $con = "1" ]; then
    if [ ! -d "$analysisdir" ]; then
	mkdir -p $analysisdir
    else
	echo "ERROR: analysis directory already existes"
	exit 1
    fi
fi

echo "starting collectl"
collectl -s cmd -f $analysisdir${processName}collectl &
CollectlPid=$!

echo "starting recording cpuMHz"
#cpuMHz.sh $analysisdir &
#cpuMHzPid=$!

echo "started monerating"

echo "starting isis"

/usr/bin/time -v docker run --privileged=true --cpuset-cpus="$processor" -v /illumina:/illumina -v /mnt/raid:/mnt/raid -v /mnt/nfs/Scritps:/mnt/nfs/Scritps ouruser/centos:plusgnuplot /illumina/development/Isis/2.6.53.8/Isis -r $workingdir -a ${analysisdir} > ${analysisdir}${processName}.stdout


echo a"Isis exit status: $?"
echo "killing monerating"
kill $CollectlPid
#kill $cpuMHzPid

collectl -p $analysisdir${processName}collectl* -P -f ${analysisdir}collectlplot

gzip -d ${analysisdir}collectlplot*

mv ${analysisdir}collectlplot* ${analysisdir}${processName}.dat

collectlTimeDateConverter.py "$processName" "$analysisdir" "$analysisdir"

cp ${analysisdir}${processName}.tsv ${analysisdir}collectl.tsv

gnuisaac.sh ${analysisdir}

rm ${analysisdir}collectl.tsv

exit
