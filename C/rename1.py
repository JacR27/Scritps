#!/usr/bin/env python
import os

def main():
    extention1 = "bcl.4bins.remapped.gz"
    extention2 = "bcl.gz"
    for lane in range(3,4):
        for folder in range(1,203):
            for serface in range(1,3):
                for swath in range(1,4):
                    for tile in range(1,9):
                        path = "./C{:d}.1/s_{:d}_{:d}{:d}{:02d}.".format(folder,lane, serface, swath,tile)
                        if os.path.isfile(path+extention1):
                            if (os.path.isfile(path+extention2))==0:
                                os.rename(path+extention1,path+extention2)
                            else:
                                print(path+extention2+" is already a file")
                                return 1
                        else:
                            print(path+extention1+" is not a file")
                            return 1

main()
                                    
