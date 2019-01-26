#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that filters an interactome file. Only pair of proteins localized
in the same subcellular compartment remains. Need the result of 
ids_location.py
Input: python3 location_filter.py <interactome_file> <Location_file1>
                <location_file2> [outname]
'''

from sys import argv,exit
import subprocess as sbp
import os

def parse_locations(loc_files):
    '''Return a dictionary {id:location} from location files 
    Id  Location
    '''
    dic={}
    for loc_file in loc_files:
        with open(loc_file) as filex:
            for line in filex:
                dic[line.split()[0]]=line.split()[1]
    return dic
    
def filt_inter_location(inter_file,loc_dic,outname):
    '''Writes a file filtering that both proteins have same location
    '''
    f=open(outname,'w')
    with open(inter_file) as filex:
        for line in filex:
            line=line.split()
            if line[0] in loc_dic and line[1] in loc_dic:
                if loc_dic[line[0]]==loc_dic[line[1]]:
                    f.write('{}\t{}\t{}\n'.format(line[0],line[1],line[2]))
    f.close()

if __name__=="__main__":
    inter_file=argv[1]
    loc_file1=argv[2]
    loc_file2=argv[3]
    if len(argv)>4:
        outname=argv[4]
    else:
        outname=inter_file.replace('.tsv','_locfilt.tsv').split('/')[-1]
    
    loc_dic=parse_locations([loc_file1,loc_file2])
    filt_inter_location(inter_file,loc_dic,outname)
