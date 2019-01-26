#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that compares two interactomes in tsv format. Pair that are not 
contained in one of them will be set to score 0
Input: python3 compare_interactomes.py <dyer_int> <real_int> [outname]
'''

from sys import argv,exit
import subprocess as sbp
import os

def parse_int(int_file):
    '''Return a dictionary with {pair:score}
    
    int_file: string, filename of interactome
    '''
    int_dic={}
    with open(int_file) as filex:
        for line in filex:
            pair=line.split()[:2]
            pair.sort()
            descr=pair[0]+'|'+pair[1]
            descr=descr.strip()
            
            int_dic[descr]=line.split()[2]
    return int_dic

def write_comparation_csv(real_dic,dyer_dic,outname):
    '''Write a csv with pair,real_score,dyer_score.
    For pairs not in the real interactome, score = 0
    '''
    f=open(outname,'w')
    for pair in dyer_dic:
        dyer_score=dyer_dic[pair]
        if pair in real_dic:
            real_score=real_dic[pair]
            f.write('{},{},{}\n'.format(pair,real_score,dyer_score))
    f.close()


if __name__=="__main__":
    dyer_int=argv[1]
    real_int=argv[2]
    if len(argv)>3:
        outname=argv[3]
    else:
        outname='pr_comp.csv'
    
    dyer_dic=parse_int(dyer_int)
    real_dic=parse_int(real_int)
    
    write_comparation_csv(real_dic,dyer_dic,outname)
    
