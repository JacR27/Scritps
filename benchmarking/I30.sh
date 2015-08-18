#!/bin/bash

export LD_LIBRARY_PATH=/home/sbsuser/local/zlibCloudFlare-install/lib
/home/sbsuser/local/iSAAC-gcc-4.7.4-install/bin/isaac-align -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format bcl-gz -m 122 --tiles s_4

exit
