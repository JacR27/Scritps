import numpy as np
import os.path
import time
import argparse
import multiprocessing
import os
import yaml
from pprint import pprint as pp 

def getoptions():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-r', required=True, action='store', dest='reference',
                            help='path to fasta reference')
    return arg_parser.parse_args()

def readfilter(filename):
    with open(filename) as fh:
        dt = np.dtype('<i4')
        np.fromfile(fh,dtype='B',count=8)
        nClusters = np.fromfile(fh,dtype=dt,count=1)
        print(nClusters)
        mask = np.fromfile(fh,dtype='B')
    return mask, nClusters
    
def readbytes(filename):
    with open(filename) as fh:
        dt = np.dtype('<i4')
        nClusters = np.fromfile(fh,dtype=dt,count=1)
        print(nClusters)
        basecalls = np.fromfile(fh,dtype='B')
    return basecalls, nClusters

def getbc(basecalls):
    bases = np.mod(basecalls, 4)
    qualities = np.divide(basecalls, 4)
    bases =bases[0]
    qualities = qualities[0]
    return qualities,bases

def read_lane_bci(filename):
    with open(filename) as fh:
        dt = np.dtype('<i4')
        tileIndexes = np.fromfile(fh,dtype=dt)
        ntiles = tileIndexes.size/2
        print(ntiles)
    for i in range(0,ntiles*2,2):
        print('Tile{}: {} clusters'.format(tileIndexes[i],tileIndexes[i+1]))

def translateQualities():
    quality_Map = [0,17,30]
    #digitize
    bins = np.array([0,1,17,30,100])
        
def printdir(folder):
    return [f for f in os.listdir(folder) if os.path.isfile(f) and os.path.splitext(f)[1]=='.bcl']
    

def bclstats(filename):
    qualities, bases = getbc(readbytes(filename))
    nbases = len(bases)
    baseRange = np.unique(bases)
    baseHist = np.bincount(bases)
    nqualities = len(qualities)
    qualityRange = np.unique(qualities)
    qualityHist = np.bincount(qualities)
    qualityHist = qualityHist[qualityRange]
    print(qualityRange, qualityHist)
    return baseRange, baseHist, qualityRange, qualityHist

def lane_bcl_stats():
    bcls = [bcl for bcl in printdir('./') if bcl.endswith('.bcl')]
    nbcls = len(bcls)
    bclmeta={}
    for bcl in bcls:
        baseRange, baseHist, qualityRange, qualityHist = bclstats(bcl)
        bclmeta[bcl] = [baseRange.tolist(), baseHist.tolist(), qualityRange.tolist(),qualityHist.tolist()]
    pp(bclmeta)
    with open('bclstats.yaml', 'w') as fh:
        yaml.dump(bclmeta,fh, default_flow_style=False)
    

def countfilter(filename):
    mask, nclusters= readfilter(filename)
    npf = np.sum(mask)
    ppf = (npf/nclusters[0])*100
    print('number of cluster: {}'.format(nclusters))
    print('number passing filter: {} ({}%)'.format(npf,ppf))

def main():
    #lane_bcl_stats()    
    countfilter('s_1.filter')

    #read_lane_bci('s_1.bci')
    #printbc(readbytes('somefile'))
    
if __name__ == '__main__':
    start_time = time.time()
    main()
    print((time.time()-start_time)/60)


        
