#!/usr/bin/env python
import os
import argparse


def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-l', required=True, type=int, action='store', dest='lane',
                        help='lane number')
    arg_parser.add_argument('-c', required=True, type=int, action='store', dest='first_cycle',
                        help='cycle start')
    arg_parser.add_argument('-i', required=True, action='store', dest='input_extention',
                        help='input extention')
    arg_parser.add_argument('-o', required=True, action='store', dest='output_extention',
                        help='output extention')
    
    arg_parser_options = arg_parser.parse_args()

    extention1 = arg_parser_options.input_extention
    extention2 = arg_parser_options.output_extention
    for lane in range(arg_parser_options.lane,arg_parser_options.lane+1):
        for folder in range(arg_parser_options.first_cycle,303):
            for serface in range(1,3):
                for swath in range(1,3):
                    for tile in range(1,25):
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
                                    
