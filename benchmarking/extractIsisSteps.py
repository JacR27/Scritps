# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 17:29:09 2015

@author: jrayner
"""
import sys

def main():
    NUMBER_OF_COLLUMN_SPLITS=4
    LAST_COLLUMN_AFTER_SPLIT=4
    
    for line in sys.stdin.readlines():
        if line.strip() != "":
            sys.stdout.write("".join([i.replace(" ","_") for i in (line.strip().split(None, NUMBER_OF_COLLUMN_SPLITS)[LAST_COLLUMN_AFTER_SPLIT])],))
            sys.stdout.write("\n")
    
main()
