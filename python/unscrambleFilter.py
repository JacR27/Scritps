#!/bin/env python
import sys

def sortTags(recored):
    tagOrder = {"SM":1, "AS":2, "RG":3, "NM":4, "BC":5, "OC":6, "SA":7}
    record = recored.strip().split()
    return "\t".join(record[:11]+sorted(record[11:], key=lambda tag: tagOrder[tag[:2]]))

def printHeader(header):
    for i in sorted(header,key=lambda line: line[:3]):
        print(i)
            
def processRecords():
    for line in sys.stdin:
        print(sortTags(line))

def processHeaderAndFirstRecord():
    header = []
    for line in sys.stdin:
        if line.startswith("@"):
            if "scramble" not in line:
                header.append(line.strip())
        else:
            printHeader(header)
            print(sortTags(line))
            return

def main():
    processHeaderAndFirstRecord()
    processRecords()

if __name__=='__main__':
    main()
