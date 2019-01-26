#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that filters a fasta file for getting only the sequences 
with descriptors listed in a txt file
Input: python3 filter_fasta.py <fasta_file.gzip> <interactions_file>
'''

from sys import argv
import gzip

def get_descriptors(interaction_file):
    '''Return a list with descriptors
    
    interaction_file: string, filename of interactors file
    '''
    descr=[]
    with open(interaction_file) as filex:
        filex.readline()
        for line in filex:
            descr+=line.split()[:2]
    descr=list(set(descr))
    descr.sort()
 
    return descr

def get_ids(filename):
    '''Return a list with protein ids
    
    filename: string, filename of ids file
    '''
    ids=[]
    
    with open(filename) as filex:
        filex.readline()
        for line in filex:
            ids+=[line.strip()]
    return ids
    
def write_filt_fasta(fasta_file, descr, inter_file):
    '''Write the filtered fasta file
    
    fasta_file: string, filename of fasta file
    descr: list, descriptors of interacting proteins
    '''
    score=inter_file.split('_')[-1].replace('.txt','')
    first=True
    copy=False
    f=open('filt_{}_'.format(score)+fasta_file.replace('.gz',''), 'w') 
    
    with gzip.open(fasta_file,'rt') as filex:
        for line in filex:
            if line.startswith('>'):
                if line[1:].strip() in descr[:]:
                    copy=True
                    descriptor=line[1:].strip()
                    if first:
                        first=False
                    else:
                        f.write('\n')
                    f.write('>'+descriptor+'\n')
                else:
                    copy=False
            else:
                if copy:
                    f.write(line.strip())
    return

if __name__=="__main__":
    fasta_file=argv[1]
    inter_file=argv[2]
    
    descr=get_descriptors(inter_file)
    
    write_filt_fasta(fasta_file,descr,inter_file)
    
