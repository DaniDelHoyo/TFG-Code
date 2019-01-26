#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that filters a GO list in a file for only those terms from a 
specific branch

Input: python3 filter_go_type.py <go_list_file> <go_branches_file> 
        branches_to_mantain
'''

from sys import argv
import subprocess as sbp

def parse_gobranches(b_file):
    '''Return a dictionary {GO:branch}
    '''
    dic={}
    with open(b_file) as filex:
        for line in filex:
            line=line.split()
            dic[line[0]]=line[1]
    return dic
    
def filter_go_file(go_file,b_dic,branches,outname):
    '''Write a file with the go file filtered
    '''
    with open(outname,'w') as f:
        with open(go_file) as filex:
            for line in filex:
                if b_dic[line.split()[0]] in branches:
                    f.write(line)
    

if __name__=="__main__":
    go_file=argv[1]
    go_branches_file=argv[2]
    if len(argv)>3:
        branches,sorti=[],''
        for i in range(3,len(argv)):
            branches+=[argv[i]]
            sorti+=argv[i].split('_')[0][0]+argv[i].split('_')[1][0]+'_'
    else:
        branches=['biological_process']
        sorti='bp_'
    outname=sorti+go_file.split('/')[-1]
    
    b_dic=parse_gobranches(go_branches_file)
    filter_go_file(go_file,b_dic,branches,outname)
    
    
