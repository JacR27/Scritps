import os

def main():
    extention1 = "bcl.gz"
    extention2 = "bcl.filtered.rm_11_34_42.gz"
    extention3 = "bcl.gz.rrr_42.gz"
    extention4 = "bcl.gz"
    for lane in range(5,6):
        for folder in range(1,303):
            for serface in range(1,3):
                for swath in range(1,3):
                    for tile in range(1,25):
                        path = "./L00{:d}/C{:d}.1/s_{:d}_{:d}{:d}{:02d}.".format(lane,folder,lane, serface, swath,tile)
                        if os.path.isfile(path+extention1):
                            if (os.path.isfile(path+extention2))==0:
                                if os.path.isfile(path+extention3):
                                    os.rename(path+extention1,path+extention2)
                                    #print(path)
                                    if os.path.isfile(path+extention3):
                                        if (os.path.isfile(path+extention4))==0:
                                            #print(path)
                                            os.rename(path+extention3,path+extention4)
                                        else:
                                            print(path+extention4+" is sudenly a file")
                                            return 1
                                    else:
                                        print(path+extention3+" is not a file anymore")
                                        return 1
                                else:
                                    print(path+extention3+" is not a file")
                                    return 1
                            else:
                                print(path+extention2+" is already a file")
                                return 1
                        else:
                            print(path+extention1+" is not a file")
                            return 1

main()
                                    
