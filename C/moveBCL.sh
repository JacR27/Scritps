#!/bin/bash

mkdir -p ${1}

rsync -a ./*.filter ${1}

for i in `seq 1 302`
do mkdir ${1}/C${i}.1
done

for i in `seq 1 302`
do mv ./C${i}.1/*${2}* ${1}/C${i}.1
done

exit
