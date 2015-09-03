#!/bin/bash 

name=$1

paulparsef.sh <$name.stdout >all$name

bamscript.sh <$name.stdout >bam$name

matchfindscript.sh <$name.stdout >match$name

