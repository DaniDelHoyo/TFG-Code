#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that run Xiaopan method to predict probability of interaction 
between proteins in two sets
Input: python3 run_xiaopan.py <domain_annot1> <motif_annot1|0>
                              <domain_annot2> <motif_annot2|0>
                              <dmdi_file> [outname]
'''

from sys import argv,exit
import subprocess as sbp
import os

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

def parse_dmdi(dmdi_file):
    '''Return a dict with {domo1|domo2:score}
    '''
    domo_dic={}
    with open(dmdi_file) as filex:
        for line in filex:
            line=line.split()
            pair=[line[0],line[1]]
            pair.sort()
            descr=pair[0]+'|'+pair[1]
            domo_dic[descr]=float(line[2])
    return domo_dic

if __name__=="__main__":
    dom1_file=argv[1]
    mot1_file=argv[2]
    dom2_file=argv[3]
    mot2_file=argv[4]
    dmdi_file=argv[5]
    if len(argv)>6:
        outname=argv[6]
    else:
        outname='pr_xiaopan.tsv'
        
    #Build annotations dictionary of host and pathogen
    anot1_dic=parse_domains(dom1_file)
    if mot1_file!='0':
        anot1_dic=parse_motifs(mot1_file,anot1_dic)
    
    anot2_dic=parse_domains(dom2_file)
    if mot2_file!='0':
        anot2_dic=parse_motifs(mot2_file,anot2_dic)
    
    #Build a list with proteins
    prots1_list=list(anot1_dic.keys())
    prots2_list=list(anot2_dic.keys())
    
    #Build DMDI scores dictionary
    domo_dic=parse_dmdi(dmdi_file)
    
    
    #################################

    pairs_done=[]
    f=open(outname,'w')
    for prot1 in anot1_dic:
        for prot2 in anot2_dic:
            pp=[prot1,prot2]
            pp.sort()
            prot_descr=pp[0]+'|'+pp[1]
#            if not prot_descr in pairs_done and prot1!=prot2:
            pairs_done+=[prot_descr]
            inter_prob=1
            for d in anot1_dic[prot1]:
                for e in anot2_dic[prot2]:
                    pair=[d,e]
                    pair.sort()
                    pair_descr=pair[0]+'|'+pair[1]
                    if pair_descr in domo_dic:
                        inter_prob*=(1-domo_dic[pair_descr])
            final_prob=1-inter_prob
            f.write('{}\t{}\t{}\n'\
            .format(prot1,prot2,str(final_prob)))
    f.close()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #
    
    
    
    
    
    
    
    
    
