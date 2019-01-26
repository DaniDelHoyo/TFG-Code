#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that merges the host interactome to the predicted host-pathogen
Input: python3 make_pathointeractome.py <host_interactome> 
                <pred_interactions> 
'''

from sys import argv,exit
import subprocess as sbp

if __name__=="__main__":
    host_file=argv[1]
    pred_files=[]
    for i in range(2,len(argv)):
        pred_files+=[argv[i]]
    if len(pred_files)==1:
        outname='inter_'+pred_files[0]
    else:
        outname=input('Outname?')
    cmd=''
    for i in pred_files:
        cmd+=' '+i
    final_cmd='cat {}{} > {}'.format(host_file,cmd,outname)
    sbp.check_call(final_cmd,shell=True)
