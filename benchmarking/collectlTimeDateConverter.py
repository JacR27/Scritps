#!/usr/bin/env python

print("hello")
from datetime import datetime
from time import mktime
import sys

def main():
    LENGTH_JUNK_FROM_START_OF_TITLE = 6
    LENGTH_OF_DATE_AND_TIME = 17
    EXTENTION_OF_INPUT_FILE = ".dat"
    EXTENTION_OF_OUTPUT_FILE = ".tsv"

    name = sys.argv[1]
    inputDirectory = sys.argv[3]
    outputDirectory = sys.argv[2]

    inputTSV = open(inputDirectory + name + EXTENTION_OF_INPUT_FILE,'r')
    outputTSV = open(outputDirectory + name + EXTENTION_OF_OUTPUT_FILE,'w')

    data = []
    for line in inputTSV: 
        if line[0]=="#":                                                  
            titleRow = line[LENGTH_JUNK_FROM_START_OF_TITLE:]
        else:
            data.append(line.strip()) 
   
    outputTSV.write(titleRow)                                            

    startTime = datetimestr2seconds(data[0][0:17])                        
    for i in data:
        time = datetimestr2seconds(i[0:LENGTH_OF_DATE_AND_TIME])-startTime
        outputTSV.write(str(time) + i[LENGTH_OF_DATE_AND_TIME:] + "\n") #write date to outputfile
    

    inputTSV.close()
    outputTSV.close()
    
def datetimestr2seconds(timeString):
    dt = datetime.strptime(timeString,"%Y%m%d %H:%M:%S") # convert time string into datatime object
    dt = int(mktime(dt.timetuple())) # convert datetime object to seconds

    return dt

main()
