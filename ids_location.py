#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that creates a file with LOCALIZER outnput that classifies the
proteins in compartments. To use by location_filter.py

Input: python3 ids_location.py <Localizer_results_folder> cytoplasm(T|[F])
'''

from sys import argv,exit
import subprocess as sbp
import os

def all_prots(res_file):
    '''Return a list with all proteins identifiers
    '''
    prots=[]
    with open(res_file) as filex:
        for line in filex:
            if not line.startswith('#'):
                if len(line.split())>0 and line.split()[0]!='Identifier':
                    prots+=[line.split()[0]]
    return prots

def write_id_location(loc_file,outname,class_prots):
    loca_file=loc_file.split('/')[-1]
    if loca_file.startswith('chloroplast'):
        loc='ch'
    elif loca_file.startswith('mitochondria'):
        loc='mt'
    elif loca_file.startswith('nucleus'):
        loc='nc'
    else:
        return class_prots
    
    f=open(outname,'a')
    with open(loc_file) as filex:
        for line in filex:
            if line.startswith(">"):
                idi=line.split()[0].replace(">","")
                class_prots+=[idi]
                f.write('{}\t{}\n'.format(idi,loc))
                
    f.close()
    return class_prots

def add_cyto(prots,class_prots,outname):
    '''Add cytoplasm location for not classified proteins
    '''
    
    with open(outname,'a') as f:
        for p in prots:
            if not p in class_prots:
                f.write('{}\tct\n'.format(p))
    

if __name__=="__main__":
    locali_folder=argv[1]
    if len(argv)>2:
        cyto=argv[2]
        outname=locali_folder+'ids_location_cyto.txt'
    else:
        cyto=False
        outname=locali_folder+'ids_location.txt'
        
    
        
    locali_files=os.listdir(locali_folder)
    
    with open(outname,'w'):
        pass
    class_prots=[]
    for resu in locali_files:
        class_prots=write_id_location(locali_folder+resu,outname,class_prots)
    
    if cyto in ['t','T','true','True']:
        prots=all_prots(locali_folder+'Results.txt')
        add_cyto(prots,class_prots,outname)
    
    
