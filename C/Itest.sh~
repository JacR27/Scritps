#!/bin/bash

/illumina/development/iSAAC/latest/bin/isaac-align --bam-gzip-level 6 --lane-number-max 8 --scatter-repeats 1 --reference-genome /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml --cleanup-intermediary 1 --memory-limit 123 --base-quality-cutoff 15 --stats-image-format none --gap-scoring bwa --variable-read-length yes --ignore-missing-bcls 1 --ignore-missing-filters 1 --split-gap-length 10000 --per-tile-tls 1 --seed-length 16 --description 150509_HSX213_0053_L005_FILTERED_NA12878 --clip-overlapping 1 -b RunInfo.xml --barcode-mismatches 1 --use-bases-mask y*,y* --clip-semialigned no --anchor-mate yes --qscore-bin no --tiles s_5 --base-calls-format bcl-gz 2>stdout.txt

exit
