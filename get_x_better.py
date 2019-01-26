#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that takes best interactions that have the scores for dyer and 
xiaopan over a threshold
Input: python3 get_better_comp.py <comp_file> <number_of_better> [<limit>] 
 '''

from sys import argv,exit
import subprocess as sbp

if __name__=="__main__":
    comp_file=argv[1]
    x_better=int(argv[2])
    if len(argv)>3:
        limit=float(argv[3])
    else:
        limit=0
        
    outname='{}_{}_{}'.\
    format('better',x_better,comp_file.split('/')[-1])
    
    sbp.check_call("awk -F ',' '{print $1,$2+$3}' "+\
    comp_file+" > pr.tmp",shell=True)
    sbp.check_call('sort -k 2 -r pr.tmp > pr2.tmp',shell=True)
    
    f=open(outname.replace('.csv','.tsv'),'w')
    with open('pr2.tmp') as filex:
        i=1
        for line in filex:
            line=line.split()
                
            pair=line[0].split('|')
            score=line[1].strip()
            if float(score)<limit:
                break
            f.write('{}\t{}\t{}\n'.format(pair[0],pair[1],score))
            i+=1
            if i>x_better:
                break
    f.close()
    sbp.check_call('rm pr.tmp',shell=True)
    sbp.check_call('rm pr2.tmp',shell=True)
