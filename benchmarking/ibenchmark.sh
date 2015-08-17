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

nohup /usr/bin/time -v $workflow.sh > ${analysisdir}${processName}.stdout

echo killing collectl

kill $CollectlPid

collectl -p $analysisdir${processName}* -P -f ${analysisdir}plot*

gzip -d ${analysisdir}plot*

mv ${analysisdir}plot* ${analysisdir}${processName}.dat

cp $validationFile ${resultsDir}${processName}.val

grep -i "error" ${analysisdir}${processName}.stdout > ${resultsDir}${processName}.error

mv $workingdir/Aligned ${processName}Aligned

echo $processName >> ${resultsDir}runtimes

grep "Elapsed (wall clock) time" ${analysisdir}${processName}.stdout >> ${resultsDir}runtimes

rm ${workingdir}/Temp/{bin-*,gnuplot-*}

scriptsDir=/home/sbsuser/Scritps/benchmarking/

grep "Elapsed time" ${analysisdir}${processName}.stdout | python ${scriptsDir}extractIsisSteps.py > ${analysisdir}$processName

python ${scriptsDir}collectlTimeDateConverter.py "$processName" "$analysisdir" "$analysisdir"

export GNUTERM=dumb
/usr/bin/gnuplot <<\EOF

reset
clear
workflow = ""
hostname = ""
file="collectl.tsv"
set terminal png size 1500,500
set output "CPU".hostname.workflow.".png"

set xdata time
set timefmt "%s"
set format x "%H:%M"
set xlabel "time"
set ylabel "Percent CPU"
#set yrange [0:110]
set title "CPU usage"." ".hostname." ".workflow
set style data line
set key outside
#set xtic (1000)

plot file using (column("Time")):(column("[CPU]Wait%")):xtic(60) title "[CPU]Wait%", "" using (column("Time")):(column("[CPU]Sys%")) title "[CPU]Sys%", "" using (column("Time")):(column("[CPU]Nice%")) title "[CPU]Nice%"# "" using (column("Time")):(column("[CPU]User%")) title "[CPU]User%"

plot file using (column("Time")):(column("[CPU]User%")) title "[CPU]User%"

replot

set output "memory".hostname.workflow.".png"
set xdata time
set timefmt "%s"
set format x "%H:%M"
set xlabel "time"
set ylabel "GB"
set title "Memory usage"." ".hostname." ".workflow
set style data line

plot file using (column("Time")):(column("[MEM]Commit")) title "[MEM]Commit", "" using (column("Time")):(column("[MEM]Cached")) title "[MEM]Cashed"

set output "IO".hostname.workflow.".png"
set xdata time
set timefmt "%s"
set format x "%H:%M"
set xlabel "time"
set ylabel "MB/s
set title "Disk IO"." ".hostname." ".workflow
set style data line
plot file using (column("Time")):(column("[DSK]ReadKBTot")) title "[DSK]ReadKBTot%", "" using (column("Time")):(column("[DSK]WriteKBTot")) title "[DSK]WriteKBTot%"
set output
#
EOF

rm collectl.tsv

exit
