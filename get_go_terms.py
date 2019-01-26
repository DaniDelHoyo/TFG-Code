#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Student number: 960824-177-040
Script that extract go terms present in proteins in a interaction file (column 0)
Input is:
    python3 get_go_terms.py <uni_interaction_file> <gaf_file> outname
'''

from sys import argv
import subprocess as sbp

def parse_gofun(gofile):
    '''Return a dic with {GOid:function}
    '''
    dic={}
    with open(gofile) as filex:
        for i in range(10):
            filex.readline()
        for line in filex:
            line=line.split('\t')
            dic[line[0]]=line[1].replace(' ','_')
    return dic
    
def parse_gaf(gaf):
    '''Returns a dictionary with {UniprotID:[GOs]}
    '''
    dic={}
    with open(gaf) as filex:
        for line in filex:
            line=line.split('\t')
            go_an=line[4]
            gen_id=line[1]
            if gen_id in dic:
                dic[gen_id]+=[go_an]
            else:
                dic[gen_id]=[go_an]
    return dic

def parse_inter(inter_file,col=0):
    '''Returns a list with proteins present in a interactome
    '''
    prots=[]
    with open(inter_file) as filex:
        for line in filex:
            line=line.split()
            if col==None:
                prots+=[line[0],line[1]]
            else:
                prots+=[line[col]]
    prots=list(set(prots))
    return prots

def filter_dic(go_dic,inter_prots):
    '''Returns a go dictionary only with proteins in inter_prots
    '''
    dic={}
    for prot in inter_prots:
        if prot in go_dic:
            dic[prot]=go_dic[prot]
    return dic

if __name__=="__main__":
    inter_file=argv[1]
    gaf_file=argv[2]
    outname=argv[3]
    
    host_prots=parse_inter(inter_file)
    gaf_dic=parse_gaf(gaf_file)
    go_func= parse_gofun('GO_functions')
    
    good_gos={}
    for p in host_prots:
        try:
            gos=gaf_dic[p]
            for go in gos:
                if go in good_gos:
                    good_gos[go][1]+=1
                else:
                    good_gos[go]=[go_func[go],1]
        except:
            pass
                
    with open(outname,'w') as f:
        for go in good_gos:
            f.write('{}\t{}\t{}\n'\
            .format(go,good_gos[go][0],good_gos[go][1]))
                
                
                
                
                
                
