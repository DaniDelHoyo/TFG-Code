#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that run netsurfp only for those proteins that are not present 
already in a netsurfp output and then merge the files
Input: python3 redo_netsurf.py <prev_nets_file> <fasta_file> <Outname>
'''

from sys import argv,exit
import subprocess as sbp
import os

def build_fasta_dic(fasta_file):
    '''Return a dictionary with {descriptor:sequence} of a fasta file
    '''
    fasta_dic={}
    with open(fasta_file) as filex:
        for line in filex:
            if line.startswith('>'):
                descr=line[1:].strip()
                fasta_dic[descr]=''
                
            else:
                fasta_dic[descr]+=line.strip()
    return fasta_dic

def write_new_fasta(fasta_dic,prev_net):
    '''Write a fasta file only with sequences not present in prev_net
    '''
    out=prev_net+'tmp.fa'
    f=open(out,'w')
    cmd='grep -o "4081.[^ ]*" {} | uniq'.format(prev_net)
    xd=str(sbp.check_output(cmd,shell=True))
    for descr in fasta_dic:
        if not descr in xd:
            f.write('>{}\n{}\n'.format(descr,fasta_dic[descr]))
    f.close()
    return out

def run_NetSurfP(new_fasta,outname):
    '''Run netsurfp with the a fasta file
    '''
    print('Running NetSurfP')
    cmd='./netsurfp -i {} -o {}'.format(new_fasta,outname)
    sbp.check_call(cmd,shell=True)
    
def merge_nets(prev_net,new_net):
    '''Merge the new and the previous nets files
    '''
    f=open(new_net,'a')
    with open(prev_net) as filex:
        for i in range(7):
            filex.readline()
        for line in filex:
            f.write(line)
    f.close()

if __name__=="__main__":
    prev_net=argv[1]
    fasta_file=argv[2]
    outname=argv[3]
    
    fasta_dic=build_fasta_dic(fasta_file)
    
    new_fasta=write_new_fasta(fasta_dic,prev_net)
    
    run_NetSurfP(new_fasta,outname)
    
    merge_nets(prev_net,outname)
    #Deleting temporary fasta')
    sbp.check_call('rm '+new_fasta,shell=True)
