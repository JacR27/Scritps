#!/bin/bash
PATH=$PATH:/mnt/nfs/Scritps/benchmarking
export PATH

usage()
{
cat <<EOF
usage:

ISIS: $ISIS

EOF
}

system=`hostname`
con=1
#get option
if [ $# -eq 0 ]; then usage; exit 1; fi
while getopts "r:w:p:d:s:m:n:c" argument; do
    case $argument in
	d ) workingdir=$OPTARG;;
        r ) run=$OPTARG;;
	w ) workflow=$OPTARG;;
	p ) processor=$OPTARG;;
	c ) con=0;;
	s ) system=$OPTARG;;
	m ) mem=$OPTARG;;
	n ) cpumem=$OPTARG;;
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
    echo "ERROR: processor (-p) required!"
    exit 1
fi
if [[ -z $mem ]]; then
    echo "ERROR: mem (-m) required!"
    exit 1
fi

if [[ -z $cpumem ]]; then
    echo "ERROR: cpumem (-n) required!"
    exit 1
else
    optmemset="--cpuset-mems=${cpumem}"
fi

if [[ -z $workingdir ]]; then
    echo "ERROR: workingdir (-d) required!"
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

/usr/bin/time -v docker run --privileged=true --cpuset-cpus="$processor" -m $mem -t $optmemset -v /illumina:/illumina -v /mnt/raid:/mnt/raid -v /mnt/ssd:/mnt/ssd -v /mnt/nfs:/mnt/nfs -v /mnt/nfs1:/mnt/nfs1 -v /mnt/nfs/Scritps:/mnt/nfs/Scritps ouruser/cenots:plusgnuplot /illumina/development/Isis/2.6.53.15/Isis -r $workingdir -a ${analysisdir}> ${analysisdir}${processName}.stdout

echo a"Isis exit status: $?"
echo "killing monerating"
kill $CollectlPid
#kill $cpuMHzPid

collectl -p $analysisdir${processName}collectl* -P -f ${analysisdir}collectlplot

gzip -d ${analysisdir}collectlplot*

mv ${analysisdir}collectlplot* ${analysisdir}${processName}.dat

collectlTimeDateConverter.py "$processName" "$analysisdir" "$analysisdir"

grep "Elapsed time" ${analysisdir}WorkflowLog.txt | extractIsisSteps.py > ${analysisdir}${processName}.stp

cp ${analysisdir}${processName}.tsv ${analysisdir}collectl.tsv

gnuisaac.sh ${analysisdir}

rm ${analysisdir}collectl.tsv

exit
