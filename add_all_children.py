#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Student number: 960824-177-040
Script that add all GO children terms to a file with GO terms
Input is:
    python3 add_all_children.py <good_gos> <children_go_file>
'''

from sys import argv
import subprocess as sbp

def take_gos(go_file,col=0):
    '''Return a list with elements in column
    '''
    listi=[]
    with open(go_file) as filex:
        for line in filex:
            line=line.split()
            listi+=[line[0]]
    return listi

def parse_chil_file(chil_file):
    '''Returns a dictionary {parent:[children]}
    '''
    dic={}
    with open(chil_file) as filex:
        for line in filex:
            line=line.split()
            dic[line[0]]=line[1:]
    return dic

def recursive_children(go_list,chil_dic):
    '''Return a list with all children of all go terms in a list
    '''
    children=[]
    for go in go_list:
        if go in chil_dic:
            children+=chil_dic[go]
        else:
            pass
            
    if len(children)>0:
        go_list+=recursive_children(children,chil_dic)
    
    return go_list

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

def write_all_gos(all_gos,outname):
    '''Write a file with original GO terms plus all their children
    '''
    #Parsing GO term with GO function
    go_func=parse_gofun('GO_functions')
    with open(outname,'w') as f:
        for go in all_gos:
            if go in go_func:
                fun=go_func[go]
            else:
                fun='?????????'
            f.write('{}\t{}\n'.format(go,fun))

if __name__=="__main__":
    go_file=argv[1]
    chil_file=argv[2]
    outname=go_file.replace('.tsv','_+children.tsv')
    
    go_list=take_gos(go_file)
    chil_dic=parse_chil_file(chil_file)
    
    all_gos=list(set(recursive_children(go_list,chil_dic)))
    
    write_all_gos(all_gos,outname)
    
