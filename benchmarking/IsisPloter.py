#!/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("resultdir", help="resultdir")
parser.add_argument("name", help="resultdir")
args = parser.parse_args()




def main():
    readCollectlData(args.name,1)


def readCollectlData(name,Smoothing):

    MINUTES = 60
    RESULTS_DIR='args.resultdir'
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


    subprocesses = open(RESULTS_DIR + name,"r") # open file containt subprocesses
    rawSubprocess= np.array([line.strip().split() for line in subprocesses])
    subprocesses.close()
    subprocessTimes = np.array(rawSubprocess[:,0],dtype=np.float)/MINUTES
    NumberOfSubprocesses=len(subprocessTimes)
    subprocessNames = np.array(rawSubprocess[:,1])
    clumaltiveTimes = []
    for i in range(NumberOfSubprocesses):
        clumaltiveTimes.append(np.sum(subprocessTimes[0:i+1])) #generate clumaltive time for ploting ticks
    clumaltiveTimes=np.array(clumaltiveTimes)

    return [columnLables, data , subprocessTimes, subprocessNames, clumaltiveTimes]

def plotRunData(runs,columns,maxTime,showMultipulLeg,runInfo):

    fig = plt.figure()
    numberOfSubplots = len(runs)
    for subplot in range(numberOfSubplots):
        ax = fig.add_subplot(numberOfSubplots,1,subplot+1)
        ax.set_xlim(0,maxTime)
        ax.set_title(runInfo[subplot].name)
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
    fig.savefig('args.resultdir' + str(1) + '.png',bbox_inches='tight')
    fig.show()




main()
