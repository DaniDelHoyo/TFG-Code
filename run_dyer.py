#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that train the probabilities following Dyer et al 2007 method
using the annotation of a interactome by domains and motifs and
the interactome itself.
Input: python3 run_dyer.py <Trained_dyer>
        <Domains_file1> <Motifs_file1>|0 
        <Domains_file2> <Motifs_file2>|0 [<final_output>]
'''

from sys import argv,exit
import subprocess as sbp
import os

def parse_args(args_file):
    arguis=[]
    with open(args_file) as filex:
        for line in filex:
            arguis.append(line.split()[-1])
    return arguis

def parse_domains(dom_file):
    '''Return a dictionary {protein:[domains]}
    
    dom_file: string, filename of domain annotation
    '''
    anot_dic={}
    with open(dom_file) as filex:
        for line in filex:
            line=line.split()
            if line[0] in anot_dic:
                anot_dic[line[0]]+=[line[4]]
            else:
                anot_dic[line[0]]=[line[4]]
    return anot_dic

def parse_motifs(mot_file,anot_dic):
    '''Return a dictionary {protein:[domains]}
    It will complement the domains already annotated
    
    mot_file: string, filename of domain annotation
    anot_dic: dictionary, {protein:[domains]}
    '''
    with open(mot_file) as filex:
        for line in filex:
            line=line.split()
            if line[0] in anot_dic:
                anot_dic[line[0]]+=[line[3]]
            else:
                anot_dic[line[0]]=[line[3]]
    return anot_dic

def save_posteriors(post_dic,outname):
    '''Write in a tsv file the result of training the Dyer model
    Dom1|Dom2   Posterior
    '''
    with open(outname,'w') as f:
        for pair in post_dic:
            f.write(pair+'\t'+str(post_dic[pair])+'\n')
    
def open_post_dic(trained_file):
    '''Return a dictionary with {prot1|prot2:posterior_prob} result of
    training Dyer method
    
    trained_file: string, file where posterior probabilities are
    '''
    post_dic={}
    with open(trained_file) as filex:
        for line in filex:
            post_dic[line.split()[0]]=float(line.split()[1])
    return post_dic
    
if __name__=="__main__":
    #Parsing arguments
    trained_file=argv[1]
    dom_file=argv[2]
    mot_file=argv[3]
    dom_pat=argv[4]
    mot_pat=argv[5]
    if len(argv)>6:
        outname=argv[6]
    else:
        outname='interact_dyer_pr.tsv'
    
    #Build annotations dictionary of host and pathogen
    host_dic=parse_domains(dom_file)
    if mot_file!='0':
        host_dic=parse_motifs(mot_file,host_dic)
        
    
    
    pat_dic=parse_domains(dom_pat)
    if mot_pat!='0':
        pat_dic=parse_motifs(mot_pat,pat_dic)
    
    
    #Build a list with proteins
    prots_host=list(host_dic.keys())
    prots_pat=list(pat_dic.keys())
    
    print('Opening trained')
    posterior_dic=open_post_dic(trained_file)
    
    ##############################################
    #Applying the Dyer method
    ##############################################
    print('Running Dyer method')
    f=open(outname,'w')
    for host_prot in host_dic:
        for pat_prot in pat_dic:
            inter_prob=1
            for d in pat_dic[pat_prot]:
                for e in host_dic[host_prot]:
                    pair=[d,e]
                    pair.sort()
                    descr=pair[0]+'|'+pair[1]
                    if descr in posterior_dic:
                        inter_prob*=(1-posterior_dic[descr])
            final_prob=1-inter_prob
            #print(host_prot,pat_prot,final_prob)
            f.write('{}\t{}\t{}\n'.format(host_prot,pat_prot,final_prob))
    f.close()

