import subprocess as subp
import sys
import argparse
import StringIO as io
import unittest

def getoptions():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--in1', required=True, dest='bam1', action='store',  help='first bam')
    arg_parser.add_argument('--in2', dest='bam2', action='store', required=True, help='second bam')    
    return arg_parser.parse_args()


def samtoolsView(bam):
    command = ['/mnt/nfs/cramming/samtools/bin/samtools', 'view']
    with open(bam) as fh:
        job = subp.Popen(command,stdin=fh,stdout=subp.PIPE,universal_newlines=True)
    print("waiting for communicate")
    return job
    
def getTags(recored):
    return recored.strip().split()[11:]


def sortTags(recored):
    record = recored.strip().split()
    return "\t".join(record[:11]+sorted(record[11:], key=keysort))
    
def main():
    args = getoptions()
    assert args.bam1.endswith('.bam')
    stream1 = samtoolsView(args.bam1)
    stream2 = samtoolsView(args.bam2)
    for line1 in stream1.readlines():
        line2 = stream2.readline()
        if not line1 == line2:
            print(line1)
            print(line2)
            return
def filter():
    command = ['/mnt/nfs/cramming/samtools/bin/samtools', 'view']
    args = getoptions()
    assert args.bam1.endswith('.bam')
    with open(args.bam1) as fh:
        job = subp.Popen(command,stdin=fh,stdout=subp.PIPE,universal_newlines=True)
    while 1:
        print(job.stdout.readline())

def sorttags():
    line = None
    while line != "":
        line = sys.stdin.readline().strip()
        if line.startswith('@'):
            print(line)
        else:
            print(sortTags(line))

def filtertags():
    line = None
    t = getTags(sys.stdin.readline())
    t = [i[:2] for i in t]
    while line != "":
        line = sys.stdin.readline().strip()
        tags = getTags(line)
        tags = [i[:2] for i in tags]
        for i,j in enumerate(tags):
            if j != t[i]:
                t.insert(i,j)

order = {"SM":1, "AS":2, "RG":3, "NM":4, "BC":5, "OC":6, "SA":7}

def keysort(tag):
    return order[tag[:2]]
    

def header(line):
    if "scramble" not in line:
        print(line.strip())
        
def printheader(header):
    for i in header:
        print(i)

def headsort(line):
    return line[:3]

def tagssunin():
    header = []
    isfirstrecord = 1
    for line in sys.stdin:
        if line.startswith("@"):
            if "scramble" not in line:
                header.append(line.strip())
        elif isfirstrecord:
            isfirstrecord = 0
            printheader(sorted(header,key=headsort))
            print(sortTags(line))
        else:     
            print(sortTags(line))


if __name__=='__main__':
    tagssunin()
