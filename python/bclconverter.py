import sys
import struct as st
import math
import multiprocessing as mp
import gzip
import os

def filterstats(SERFACES,SWATHS,TILES):
    filters = []
    for surface in SERFACES:
        for swath in SWATHS:
            for tile in TILES:
                filename= "s_4_"+ surface + swath + tile + ".filter"
                print(filename)
                file = open(filename,"rb")
                file.read(8)
                numberClusters = st.unpack("I",file.read(4))[0]
                bytes = file.read()
                file.close()
                data = st.unpack(str(numberClusters)+"B",bytes)
                print(sum(data))
                filters.append(data)
    return filters


def bclConverter(args):    
    def readBCL(foldername, filename, fileExtention):
        file = gzip.open(foldername+filename+fileExtention,"rb")
    
        hbytes = file.read(4) # remove first 4 bytes
        nClusters = st.unpack('I',hbytes)[0]
        dbytes = file.read()
        print(foldername + filename)
        data = st.unpack(str(nClusters)+'B',dbytes)
    
        file.close()
        return nClusters, data

    def readFasterq(filename,readlen):
        file = gzip.open(filename,"rb")
        hbytes = file.read(4)
        nClusters = st.unpack('I',hbytes)[0]
        dbytes = file.read()
        data = st.unpack(str(nClusters*readlen)+'B',dbytes)
        return nClusters, data

    def filterData(data,filters):
        newClusters = sum(filters)
        filteredData = [point for i, point in enumerate(data) if filters[i]]
        return filteredData, newClusters
        
    def extractBQ(data):
        bases =  [i%4 for i in data]
        qualities= [i//4 for i in data]
        return qualities, bases
    
    def remapQualities(qualities,qualityMap):
        remapedQualities = [qualityMap[i] for i in qualities]
        return remapedQualities

    def join(array, bytesToJoin):
        lengthOfNewArray = int(math.ceil(len(array)/float(bytesToJoin)))
        pShift = int(8//bytesToJoin)
        joined = [0] * lengthOfNewArray
        maxMod = bytesToJoin -1
        for i in range(len(array)):
            newindex = int(i//bytesToJoin)
            shift = 2**((maxMod-(i%bytesToJoin))*pShift)
            joined[newindex] = joined[newindex] + array[i]*shift
        return joined
            
    def saveArray(data,header, foldername, filename,fileExtention):
        hbytes = st.pack('I',header) 
        dbytes = st.pack(str(len(data)) +'B',*data)
        outfile=gzip.open(foldername+filename+fileExtention,'w')
        outfile.write(hbytes)
        outfile.write(dbytes)
    
    def interleave(array1,array2):
        interleved = [j for i in zip(array1,array2) for j in i]
        return interleved
        
    def basicConverter(foldername,filename):    
        OrigonalClusters, data = readBCL(foldername, filename)
        qualities, bases = extractBQ(data)
        qualityMap = {0:0, 7:1, 11:1, 22:2, 27:2, 32:2, 37:3, 42:3}
        remapedQualities = remapQualities(qualities, qualityMap)
        interleved = interleave(remapedQualities,bases)
    
        joined = join(interleved,4)
        saveArray(joined,OrigonalClusters,foldername,filename,".rqb.gz")

    def basicConverter(foldername,filename):
        OrigonalClusters, data = readBCL(foldername, filename, ".fbcl.gz")
        qualities, bases = extractBQ(data)
        qualityMap = {0:0, 7:1, 11:1, 22:2, 27:2, 32:2, 37:3, 42:3}
        remapedQualities = remapQualities(qualities, qualityMap)
        interleved = interleave(remapedQualities,bases)
        joined = join(interleved,4)
        saveArray(joined,OrigonalClusters,foldername,filename,".frqb.gz")

    def filterAndConvert(foldername,filename,filters):
        OrigonalClusters, data = readBCL(foldername, filename)
        data,newClusters = filterData(data,filters)
        saveArray(data,newClusters,foldername,filename,".fbcl.gz")
        
    def demultiplex(filename,CYCLES,readnum):
        read = []
        for cycle in CYCLES:
            foldername = "./C" + cycle + ".1/"
            OrigonalClusters, data =  readBCL(foldername, filename)
            read.append(data)
        readTf = flatten2d(transpose(read))
        saveArray(readTf,OrigonalClusters,"./",filename+readnum,".fasterq.gz")

    def demultiplexAndFilter(filename,CYCLES,readnum,filters):
        read = []
        for cycle in CYCLES:
            foldername = "./C" + cycle + ".1/"
            OrigonalClusters, data =  readBCL(foldername, filename)
            read.append(data)
        readTf, newClusters = filterData(transpose(read),filters)
        readTf = flatten2d(readTf)
        saveArray(readTf,newClusters,"./",filename+readnum,".ffasterq.gz")

    def convertFastq(filename,readlen):
        OrigonalClusters, data = readFasterq(filename+".fasterq.gz",readlen)
        qualities, bases = extractBQ(data)
        qualityMap = {0:0, 7:1, 11:1, 22:2, 27:2, 32:2, 37:3, 42:3}
        remapedQualities = remapQualities(qualities, qualityMap)
        interleved = interleave(remapedQualities,bases)
        joined = join(interleved,4)
        saveArray(joined, OrigonalCllusters, "./", filename, ".rfasterq.gz")
        
    def transpose(array):
        transposedRead = list(map(list,zip(*array)))
        return transposedRead
    
    def flatten2d(array):
        return [i for sublist in array for i in sublist]

    basicConverter(*args)
    
if __name__ == "__main__":
    
    CYCLES = tuple(str(i) for i in range(1,303))
    SERFACES = ("1","2")
    SWATHS = ("1","2")
    TILES = tuple("{:02n}".format(i) for i in range(1,25))
    
    for cycle in CYCLES:
        pool = mp.Pool()
        for serface in SERFACES:
            for swath in SWATHS:
                for tile in TILES:
                    foldername = "./C" + cycle + ".1/"
                    filename ="s_4_" + serface + swath + tile 
                    #bclConverter([foldername, filename])
                    #pool.apply_async(bclConverter,args=([filename,CYCLES[0+(int(read)-1)*151:151+(int(read)-1)*152],read,filters[sst]],))
                    pool.apply_async(bclConverter,args=([foldername,filename],))
        pool.close()
        pool.join()

