#!/bin/bash

export LD_LIBRARY_PATH=/home/sbsuser/local/gcc-5.2.0-install/lib:/home/sbsuser/local/gcc-5.2.0-install/lib64
/home/sbsuser/local/iSAAC-gcc-5.2.0-O3-install/bin/isaac-align -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format bcl-gz -m 122 --tiles s_4

exit
