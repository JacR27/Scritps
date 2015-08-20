#!/bin/bash

/illumina/development/iSAAC/iSAAC-03.15.07.30/bin/isaac-align -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format fastq-gz -m 122 --tiles s_4 --base-calls-directory /mnt/141128_HSX177_0082_BH00UABBXX/Fastq

exit
