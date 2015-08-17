#!/bin/bash

processName=$1

outputDir=/home/jrayner/benchmarking/results/
inputDir=./
scriptsDir=/home/jrayner/benchmarking/scripts/

grep "Elapsed time" ${inputDir}${processName}.stdout | python ${scriptsDir}extractIsisSteps.py > $outputDir$processName

python ${scriptsDir}collectlTimeDateConverter.py "$processName" "$outputDir" "$inputDir"
echo $outputDir
echo count > $outputDir/count
