#!/bin/env python

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import numpy as np
from datetime import datetime
from time import mktime

parser = argparse.ArgumentParser()
parser.add_argument("resultdir", help="resultdir")
parser.add_argument("name", help="name")
args = parser.parse_args()




def main():
    dataToPlot =   [["[CPU]Wait%","[CPU]Nice%","[CPU]User%","[CPU]Sys%"],
                    ["[DSK]ReadKBTot","[DSK]WriteKBTot"],
                    ["[MEM]Cached","[MEM]Commit","[MEM]Tot","[MEM]Anon"]]
    
    columns=[["[CPU]Wait%","[CPU]Nice%","[CPU]User%","[CPU]Sys%"]]

    plotRunData([truncate(*readCollectlData(args.name,10000),from_="Statistics_evaluation",to="Report_generation")],columns,1,1,1)

def readCollectlData(name,Smoothing):
    MINUTES = 60
    RESULTS_DIR=args.resultdir
    InputFileExtention = ".tsv"
    TITLE_ROW_NUMBER=1

    titleRow = open(RESULTS_DIR + name + InputFileExtention,"r")
    columnLables = titleRow.readline().split() # read colmn headding in to array
    titleRow.close()

    data = (np.loadtxt(RESULTS_DIR+name+ InputFileExtention,skiprows=TITLE_ROW_NUMBER)) #read system into numpy array
    data[:,0] = data[:,0]/MINUTES
    dataTemp = np.array(data[:,1:])
    for d in range(np.size(data,axis=0)):
        data[d,1:] = np.mean(dataTemp[max([0,d-Smoothing]):d+1,:],axis=0)


    subprocesses = open(RESULTS_DIR + name + ".stp","r") # open file containt subprocesses
    rawSubprocess= np.array([line.strip().split() for line in subprocesses])
    for i in range(np.size(rawSubprocess,0)):
        rawSubprocess[i,0]=  ts2s(rawSubprocess[i,0])
    subprocesses.close()
    subprocessTimes = np.array(rawSubprocess[:,0],dtype=np.float)/MINUTES
    NumberOfSubprocesses=len(subprocessTimes)
    subprocessNames = np.array(rawSubprocess[:,1])
    clumaltiveTimes = []
    for i in range(NumberOfSubprocesses):
        clumaltiveTimes.append(np.sum(subprocessTimes[0:i+1])) #generate clumaltive time for ploting ticks
    clumaltiveTimes=np.array(clumaltiveTimes)

    return [columnLables, data , subprocessTimes, subprocessNames, clumaltiveTimes]

def ts2s(timeString):
    dt = datetime.strptime(timeString[0:8],"%H:%M:%S") # convert time string into datatime object
    dt = int(mktime(dt.timetuple())) # convert datetime object to seconds

    return dt

def plotRunData(runs,columns,maxTime,showMultipulLeg,runInfo):

    fig = plt.figure()
    numberOfSubplots = len(runs)
    for subplot in range(numberOfSubplots):
        ax = fig.add_subplot(numberOfSubplots,1,subplot+1)
        ax.set_xlim(0,3)
        ax.set_ylim(0,3)
        #ax.set_title(runInfo[subplot].name)
        run = runs[subplot]
        for collumn in columns[subplot]:
            ax.plot(run[1][:,0],run[1][:,run[0].index(collumn)],label=collumn)
        for i , subprocess in enumerate(run[3]):
            ymin, ymax= ax.get_ylim()
            ax.annotate(subprocess, xy = (run[4][i],1), xytext = (run[4][i],ymax/1.2),rotation=60, arrowprops=dict(visible=True, fill=False, width=0.0001,linestyle='dashed'), fontsize=8)
        if showMultipulLeg:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    if not showMultipulLeg:
        ax.legend(loc='center left', bbox_to_anchor=(1, 1))
    ax.set_xlabel("Time (minutes)")
    fig.tight_layout()
    fig.set_size_inches(7,6)
    fig.savefig(args.resultdir + str(1) + '.png',bbox_inches='tight')
    fig.show()

def truncate(columnLables,data,subprocessTimes,subprocessNames,clumaltiveTimes,from_="start",to="end"):
    times = data[:,0]
    lenghtTimes = len(times)
    toTime = times[-1]
    fromTime = 0
    subFromIndex = -1
    subToIndex = len(subprocessNames)
    subBaseTime = 0
    for i,name in enumerate(subprocessNames):
        if name ==  to:
            toTime = clumaltiveTimes[i]
            subToIndex = i+1
        if name == from_:
            fromTime = clumaltiveTimes[i]
            subFromIndex = i
            subBaseTime = clumaltiveTimes[subFromIndex]
    fromIndex = 0
    toIndex = lenghtTimes-1
    
    for i,time in enumerate(times):
        if time < fromTime:
            fromIndex = i
        if times[(lenghtTimes-1)-i] > toTime:
            toIndex = (lenghtTimes-1)-i
    
    baseTime= times[fromIndex]
    subprocessNames = subprocessNames[subFromIndex+1:subToIndex]
    clumaltiveTimes = clumaltiveTimes[subFromIndex+1:subToIndex] - subBaseTime
    subprocessTimes = subprocessTimes[subFromIndex+1:subToIndex]
    data = data[fromIndex:toIndex,:]
    data[:,0] = data[:,0] - baseTime
    
    return [columnLables,data,subprocessTimes,subprocessNames,clumaltiveTimes]



main()
