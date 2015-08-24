#!/bin/bash

echo 200000 >/proc/sys/vm/nr_overcommit_hugepages
echo 40 >/proc/sys/vm/nr_hugepages

LD_LIBRARY_PATH=/usr/lib64 HUGETLB_MORECORE=yes LD_PRELOAD=/usr/lib64/libhugetlbfs.so /home/sbsuser/local/iSAAC-gcc-4.7.4-install/bin/isaac-align -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format bcl-gz -m 122 --tiles s_4

echo 0 >/proc/sys/vm/nr_overcommit_hugepages
echo 0 >/proc/sys/vm/nr_hugepages

exit
