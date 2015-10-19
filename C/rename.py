import os

extention1 = "bcl.gz"
extention2 = "bcl"
#extention3 = "bcl.gz.rrr"

for lane in range(5,6):
    for folder in range(1,303):
        for serface in range(1,3):
            for swath in range(1,3):
                for tile in range(1,25):
                    path = "./L00{:d}/C{:d}.1/s_{:d}_{:d}{:d}{:02d}.".format(lane,folder,lane, serface, swath,tile)
                    os.rename(path+extention1,path+extention2)
                    #os.rename(path+extention3,path+extention1)

