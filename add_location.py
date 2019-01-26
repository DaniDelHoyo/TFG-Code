#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that add the predicted location to a interactome file
Input:
    python3 add_location.py <interaction_file> <ids_location_host>
'''

from sys import argv,exit
import subprocess as sbp

def parse_location(loc_file):
    '''Return a dictionary with {id:location}
    '''
    dic={}
    with open(loc_file) as filex:
        for line in filex:
            line=line.split()
            dic[line[0]]=line[1]
    return dic

def write_loc_interactions(int_file,loc_dic):
    '''Write a file with the interactions and their location
    '''
    outname=int_file.replace('.tsv','_loc.tsv')
    with open (outname,'w') as f:
        with open (int_file) as filex:
            for line in filex:
                if line.split()[0] in loc_dic:
                    f.write(line.strip()+'\t'+loc_dic[line.split()[0]]+'\n')

if __name__=="__main__":
    int_file=argv[1]
    loc_file=argv[2]
    
    loc_dic=parse_location(loc_file)
    write_loc_interactions(int_file,loc_dic)
