#!/usr/bin/env python3


'''
Author: Daniel Del Hoyo Gomez
Student number: 960824-177-040
Script that quit redundancies in an .actions file from STRING (gzipped)
Option for filering by score threshold and interaction type
Input is:
    python3 deredundant.py \
    <.actions_file(e.g. 4787.protein.actions.v10.5.txt.gz)> \
    [<interaction_type(e.g. binding)>]
    [score threshold(0<th<1000)]

'''

from sys import argv
import subprocess as sbp
import os
import gzip

def build_dicc(param):
    '''Return a dictionary {prot1|prot2:score}
    
    param: dictionary with keys:
    -filename: string, filename of the .actions file
    -score_threshold: integer, score threshold for an interacton 
    -interaction_type: string, use only this type of interaction
    
    If the dictionary already has {prot2|prot1:score}, 
    {prot1|prot2:score} is not added
    '''
    dicc={}
    with gzip.open(param['filename'],'rt') as filex:
        filex.readline()
        for line in filex:
            prot1=line.split()[0]
            prot2=line.split()[1]
            score=line.split()[-1]
            if param['interaction_type'] in line and\
             int(score)>=param['score_threshold']:
                desc1=prot1+'|'+prot2
                desc2=prot2+'|'+prot1
                if desc1 in dicc or desc2 in dicc:
                    pass
                else:
                    dicc[desc1]=score
    return dicc
    
def write_nr_file(param, dicc):
    '''Write the non-redundant file
    
    filename: string, filename of the interactions file
    dicc: dictionary {prot1|prot2:score}
    '''
    filename='nr_{}_{}_{}.txt'.format(param['filename'].split('.')[0],\
    param['interaction_type'],param['score_threshold'])
    with open(filename,'w') as f:
        f.write('Accession_ID_1\tAccession_ID_2\tScore\n')
        for key in dicc:
            prot1,prot2=key.split('|')
            f.write('{}\t{}\t{}\n'.format(prot1,prot2,dicc[key]))


if __name__=="__main__":
    
    parameters={}
    parameters['filename']=argv[1]
    parameters['interaction_type']='.'
    parameters['score_threshold']=150
    
    for i in range(2,len(argv)):
        try:
            parameters['score_threshold']=int(argv[i])
        except:
            parameters['interaction_type']=argv[i]
        
    
    
    dicc=build_dicc(parameters)
    write_nr_file(parameters,dicc)
