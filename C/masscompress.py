import os
import subprocess



def main():
    extention1 = "bcl.filtered.4bins.remapped.gz.packed.gz"
    extention2 = "bcl.filtered.4bins.remapped.packed.gz"
    for lane in range(5,6):
        for folder in range(1,303):
            for serface in range(1,3):
                for swath in range(1,3):
                    for tile in range(1,25):
                        path = "./L00{:d}/C{:d}.1/s_{:d}_{:d}{:d}{:02d}.".format(lane,folder,lane, serface, swath,tile)
                        subprocess.call("")

main()
                                    
