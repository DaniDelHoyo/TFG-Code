#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that creates or append a DDI database in a specific format to use
Adapt DDI database from 3did in a common format
Input: python3 build_DDI_database.py <3did_DDI_file> 
    [out_name=DDI_database.txt]
'''

from sys import argv,exit
import subprocess as sbp
import os
import re

def ddi_parser(ifile):
    '''Parse the 3did file and return a list [[pfam1,pfam2,zscore],[]]
    
    ifile: string, filename of 3did DDI 
    '''
    inters=[]
    tot_scores=[]
    with open(ifile) as filex:
        for line in filex:
            if line.startswith('#=ID'):
                line=re.sub('\.\d+@Pfam','',line)
                line=line.split()
                pair=[line[3].replace('(',''),line[4].replace(')','')]
                scores=[]
            elif line.startswith('#=3D'):
                scores.append(float(line.split()[5]))
                tot_scores.append(float(line.split()[5]))
            elif line.startswith('//'):
                pair.append(sum(scores)/len(scores))
                inters+=[pair]
    return inters,tot_scores

def write_DDI_file(inters,tot_scores,outname):
    '''Write a tsv file Pfam1, Pfam2, score/max_score
    
    inters: list of lists [[pfam1,pfam2,score],[]]
    tot_scores: list with all the scores
    '''
    maxi_score=max(tot_scores)
    #We substract the minimum value because there are negatives
    min_score=min(tot_scores)
    with open(outname,'w') as f:
        for pair in inters:
            f.write('{}\t{}\t{}\n'.\
            format(pair[0],pair[1],(pair[2]-min_score)/maxi_score))

if __name__=="__main__":
    ifile=argv[1]
    if len(argv)>2:
        outname=argv[2]
    else:
        outname='DDI_database.txt'
    inter_list,tot_scores=ddi_parser(ifile)
    write_DDI_file(inter_list,tot_scores,outname)
    
    
