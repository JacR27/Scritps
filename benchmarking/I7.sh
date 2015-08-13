#!/bin/bash

/illumina/development/iSAAC/iSAAC-03.15.07.30/bin/isaac-align -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format bcl-gz -m 122 --tiles s_4 --temp-concurrent-load 3

exit
