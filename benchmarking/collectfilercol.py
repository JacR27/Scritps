#!/bin/env python

import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', action="store_true", help="print title line and exit")
    parser.add_argument('-n', action="store", help="collumn number or name")
    args = parser.parse_args()

    lines = sys.stdin.readlines()
    firstline=-1
    for line in lines:
        if not line.startswith('#'):
            break
        firstline += 1
        
    for i in range(firstline,len(lines)):
        line=lines[i].strip().split()
        if args.t:
            for i, j in enumerate(line):
                print("{}: {}".format(i,j))
            return
        elif args.n:
            if args.n.isdigit():
                print(line[int(args.n)])
            else:
                dex = line.index(args.n)
                print(line[dex])
        else:
            print(" ".join(line))


if __name__ == "__main__":
    main()
