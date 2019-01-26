#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that creates two tab files: 
    Parent_terms: GO_id\tParent_term1\Parent_term2...
    Child_terms:  GO_id\tChild_term1\Child_term2...
Input: python3 parse_obo.py <obo_file>
'''

from sys import argv
import subprocess as sbp

def parse_obo_dic(obo):
    '''Return a dictionary {GO_id:[GO_parents]}
    '''
    dic={}
    with open(obo) as filex:
        for line in filex:
            if line.startswith('id: GO:'):
                ide=line.split(' ')[1].strip()
                dic[ide]=[]
            elif line.startswith('is_a: GO:'):
                dic[ide]+=[line.split(' ')[1].strip()]
    return dic

def get_children(go_dic):
    '''Return a dictionary {GO_id:[GO_children]}
    '''
    chil_dic={}
    for chil in go_dic:
        for par in go_dic[chil]:
            if par in chil_dic:
                chil_dic[par]+=[chil]
            else:
                chil_dic[par]=[chil]
    return chil_dic
    
def parse_type(obo):
    '''Return a dictionary with {GO:type}
    '''
    dic={}
    with open(obo) as filex:
        for line in filex:
            if line.startswith('id: '):
                go=line.split()[-1]
                dic[go]=''
            elif line.startswith('namespace:'):
                typex=line.split()[-1]
                dic[go]=typex
    return dic
    
def write_dic(dic,outname):
    '''Write a tab file
    '''
    with open(outname,'w') as f:
        for go in dic:
            f.write(go)
            for par in dic[go]:
                f.write('\t'+par)
            f.write('\n')

def write_type(dic,outname):
    '''Write a tab file
    '''
    with open(outname,'w') as f:
        for go in dic:
            f.write('{}\t{}\n'.format(go,dic[go]))
            

if __name__=="__main__":
    go_parents_dic=parse_obo_dic(argv[1])
    go_children_dic=get_children(go_parents_dic)
    
    type_dic=parse_type(argv[1])
    
    write_dic(go_parents_dic,'pr')
    write_dic(go_children_dic,'pr2')
    write_type(type_dic,'GO_branches.tsv')
    
    sbp.check_call('sort pr > parent_gos.tsv',shell=True)
    sbp.check_call('sort pr2 > children_gos.tsv',shell=True)
    sbp.check_call('rm pr pr2',shell=True)
