#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that creates or append a DMI database in a specific format to use
Adapt DMI databases from 3did and elm in a common format
Input: python3 build_DMI_database.py <3did_DMI_file> <elm_DMI_file>
        [out_name=DMI_database.txt]
'''

from sys import argv,exit
import subprocess as sbp
import os

def build_pfam_dic(pfam_file):
    '''Return a dictionary with {identifier:accession number}
    
    pfam_file: string, filename of pfam dictionary
    '''
    pfam_dic={}
    with open(pfam_file) as filex:
        for line in filex:
            acces=line.split()[0]
            ident=line.split('\t')[-2]
            pfam_dic[ident]=acces
    return (pfam_dic)

def did_parser(did_file,pfam_dic,out_name):
    '''Write in a file the DMI in tsv format from 3did
    Accession_domain    Identifier_domain   Motif_name
    '''
    f=open(out_name,'a')
    with open(did_file) as filex:
        for line in filex:
            if line.startswith("#=ID"):
                line=line.split()
                acces=pfam_dic[line[1]]
                
                f.write(acces+'\t'+line[1]+'\t'+line[3]+'\n')
    f.close()
    return out_name 

def elm_parser(elm_file,out_name):
    '''Write in a file the DMI in tsv format from ELM
    Accession_domain    Identifier_domain   Motif_name
    '''
    f=open(out_name,'a')
    with open(elm_file) as filex:
        filex.readline()
        for line in filex:
            if "PF" in line:
                line=line.replace('"','').split()
                f.write(line[1]+'\t'+line[2]+'\t'+line[0]+'\n')
    f.close()
    return out_name

if __name__=="__main__":
    #Naming files
    did_file=argv[1]
    elm_file=argv[2]
    pfam_file=argv[3]
    if len(argv)>4:
        out_namef=argv[4]
    else:
        out_namef='DMI_database.txt'
    out_name='tmp_DMI.txt'
    
    #Creating the Pfam names dictionary
    pfam_dic=build_pfam_dic(pfam_file)
    
    #Checking if output file already exists and overwrite
    if os.path.exists(out_name):
        with open(out_namef,"w") as f:
            pass
            
    #Parsing DMI files and writing database
    did_parser(did_file,pfam_dic,out_name)
    elm_parser(elm_file,out_name)
    
    #Sorting and renaming database
    sbp.check_call('sort -k 1 tmp_DMI.txt > '+out_namef,shell=True)
    sbp.check_call('rm tmp_DMI.txt',shell=True)
    
