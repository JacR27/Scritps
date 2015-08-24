#!/bin/bash

/home/sbsuser/local/iSAAC-gcc-4.7.4-install/bin/isaac-align -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format fastq-gz -m 122 --tiles s_4 --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq --clusters-at-a-time 4000000

exit
