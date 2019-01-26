#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that mixes the DDI and the DMI databases in an unic file and 
gives a arbitrary probability to DMIs.
Sort the pairs
Input: python3 mix_DMDI.py <DDI_file> <DMI_file> <DMI_score>
        [out_name=DMDI_database.tsv]
'''

from sys import argv,exit
import subprocess as sbp
import os

def parse_DDI(ddi_file,outname):
    '''Add DDI interactions to file
    '''
    f=open(outname,'w')
    with open(ddi_file) as filex:
        for line in filex:
            line=line.split()
            pair=[line[0],line[1]]
            pair.sort()
            
            f.write('{}\t{}\t{}\n'.format(pair[0],pair[1],line[2]))
    f.close()

def parse_DMI(dmi_file,dmi_score,outname):
    '''Add DMI interactions to file with the arbitrary score
    '''
    f=open(outname,'a')
    with open(dmi_file) as filex:
        for line in filex:
            line=line.split()
            pair=[line[0],line[2]]
            pair.sort()
            
            f.write('{}\t{}\t{}\n'.format(pair[0],pair[1],dmi_score))
    f.close()

if __name__=="__main__":
    ddi_file=argv[1]
    dmi_file=argv[2]
    dmi_score=argv[3]
    if len(argv)>4:
        outname=argv[4]
    else:
        outname='DMDI_{}_database.tsv'.format(dmi_score)
    
    parse_DDI(ddi_file,outname)
    if dmi_file!='0':
        parse_DMI(dmi_file,dmi_score,outname)
    
    
    
