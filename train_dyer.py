#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that train the probabilities following Dyer et al 2007 method
using the annotation of a interactome by domains and motifs and
the interactome itself. Splicing only used for spliced interactomes.
Input: python3 train_dyer.py <interactome_file> <Domains_file> 
        <Motifs_file>|0 splicing(T/F) [<final_output>]
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

def parse_int(int_file):
    '''Return a list of tuples [(prot1,prot2),()]
    
    int_file: string, name of interactome file
    '''
    pairs=[]
    with open(int_file) as filex:
        filex.readline()
        for line in filex:
            pairs.append(tuple(line.split()[0:2]))
    return pairs


def calc_Ps(anot_dic, d,e):
    '''Return:
        pd=number of proteins with domain d
        pe=number of proteins with domain e
        pde=number of proteins with both domains
        
    '''
    pd,pe,pde=0,0,0
    for prot in anot_dic:
        ese=False
        if d in anot_dic[prot]:
            pd+=1
            ese=True
        if e in anot_dic[prot]:
            pe+=1
            if ese:
                pde+=1
    return pd,pe,pde

def calc_S(host_dic,pairs,d,e):
    '''Return the number of interactions where one protein has domain d 
    and the other has domain e
    
    '''
    
    sde=0
    for pair in pairs:
        #Checking if the proteins are annotated
        if pair[0]!=pair[1] and pair[0] in host_dic and pair[1] in host_dic:
            #Checking if the domains are in the annotated proteins
            if d in host_dic[pair[0]] and e in host_dic[pair[1]]:
                sde+=1
            elif e in host_dic[pair[0]] and d in host_dic[pair[1]]:
                sde+=1
    return sde

def calc_posterior(host_dic,pairs,d,e):
    '''Return the posterior probability given Dyer formula
    '''
    sde=calc_S(host_dic,pairs,d,e)
    pd,pe,pde=calc_Ps(host_dic,d,e)
    
    if (pd*pe)-pde!=0:
        post=sde/((pd*pe)-pde)
    else:
        post=0
        
    return post

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

def fuse_splicings(host_dic):
    '''If there are splicng forms (ATG..._1/2), the gen annotation is
    added to the dictionary as the su of al domains of splicing forms
    '''
    host_prots=list(host_dic.keys())
    splis={}
    h_dic={}
    for i in host_prots:
        if '_' in i:
            gen=i.split('_')[0]
            if gen in splis:
                splis[gen]+=[i]
            else:
                splis[gen]=[i]
        else:
            h_dic[i]=host_dic[i]
            
    for gen in splis:
        h_dic[gen]=[]
        for spli in splis[gen]:
            h_dic[gen]+=host_dic[spli]
        h_dic[gen]=list(set(h_dic[gen]))
    return h_dic

if __name__=="__main__":
    int_file=argv[1]
    dom_file=argv[2]
    mot_file=argv[3]
    fuse=argv[4]
    if fuse=='T' or fuse=='True':
        fuse=True
    else:
        fuse=False
    if len(argv)>5:
        outname=argv[5]
    else:
        outname='trained_dyer_pr.tsv'
        
    #Parsing interactome
    pairs=parse_int(int_file)
        
    #Build annotations dictionary of host and pathogen
    host_dic=parse_domains(dom_file)
    if mot_file!='0':
        host_dic=parse_motifs(mot_file,host_dic)
    
    #If the the annotation was done over spliced forms that doesn't appear 
    # in the interactome file, the annotations of those spliced forms
    # are merged with the gen name
    fuse_splis=fuse
    if fuse_splis:
        host_dic=fuse_splicings(host_dic)
    
    #A list with domains and motifs of host
    dm=list(host_dic.values())
    domos=[]
    for i in dm:
        domos+=i
    host_domos=list(set(domos))

    ##############################################
    #Starting training the model 
    #############################################
    #Training
    tot=len(host_domos)^2
    print(tot,'domains to train')
    posterior_dic={}
    for d in host_domos:
        for e in host_domos:
            #Not MMIs
            if 'PF' in d or 'PF' in e:
                domis=[d,e]
                domis.sort()
                descr=domis[0]+'|'+domis[1]
                
                post=calc_posterior(host_dic,pairs,d,e)
                posterior_dic[descr]=post
        print(d,' complete')
            
    #Saving results in a file
    save_posteriors(posterior_dic,outname)
