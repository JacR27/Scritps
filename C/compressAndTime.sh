#!/bin/bash

/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered -h bcl.filtered.1hf.gz -a 1 -g 2 &> bclfiltered1hf.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered -h bcl.filtered.9rle.gz -a 9 -g 3 &> bclfiltered9rle.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered.4bins -h bcl.filtered.4bins.9hf.gz -a 9 -g 2 &> bclfiltered4bins9hf.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered.4bins -h bcl.filtered.4bins.9rle.gz -a 9 -g 3 &> bclfiltered4bins9rle.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered.4bins.remapped -h bcl.filtered.4bins.remapped.9hf.gz -a 9 -g 2 &> bclfiltered4binsremapped9hf.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered.4bins.remapped -h bcl.filtered.4bins.remapped.9rle.gz -a 9 -g 3 &> bclfiltered4binsremapped9rle.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered.4bins.remapped.packed -h bcl.filtered.4bins.remapped.packed.9hf.gz -a 9 -g 2 &> bclfiltered4binsremappedpacked9hf.stdout

#/usr/bin/time -v BCLcompresser -l 5 -c 1 -z bcl.filtered.4bins.remapped.packed -h bcl.filtered.4bins.remapped.packed.9rle.gz -a 9 -g 3 &> bclfiltered4binsremappedpacked9rle.stdout

exit
