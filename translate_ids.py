#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Student number: 960824-177-040
Script that translate string ids in a file into uniprot ids and 
creates a new file
Input is:
    python3 translate_ids.py <interactions_file> <transaltions_file>
'''

from sys import argv
import subprocess as sbp

def build_dict(file_dict):
    '''Return a dictionary {string_id:uniprot_id}
    If there are several uniprot ids, the one with bigger score is chosen
    
    file_dic: string, filename where equivalencies are
    '''
    id_dic={}
    score_dic={}
    with open (file_dict) as filex:
        for line in filex:
            uniprot=line.split()[1].split("|")[0]
            if line.split()[0]=='4787':
                string=line.split()[2][:-2]
            else:
                string=line.split()[0]+'.'+line.split()[2]
            try:
                score=float(line.split()[-1])
                if not string in id_dic:
                    id_dic[string]=uniprot
                    score_dic[string]=score
                else:
                    if score_dic[string]<score:
                        id_dic[string]=uniprot
                        score_dic[string]=score
            except:
                pass
    if extra!=None:
        with open(extra) as filex:
            for line in filex:
                line=line.split()
                id_dic[line[2]]=line[0]
    return id_dic

def write_interactions(interaction_file, dicc,force_write=True):
    '''Write a file with translated interactions
    
    interaction_file: string, filename where interactions are
    dicc: dictionary {string_id:uniprot_id}
    
    If any string_id is not recognised, it's stored in other file
    '''
    int_name=interaction_file.split('/')[-1]
    ft=open('uni_'+int_name,'w')
    #fu=open('untrans_'+int_name,'w')
    untrans=[]
    with open(interaction_file) as filex:
        for line in filex:
            #print(line)
            str_id1=line.split()[0]
            str_id2=line.split()[1]
            if len(line.split())>2:
                int_score=line.split()[2]
            else:
                int_score=1
            if not str_id1 in dicc:
                untrans+=[str_id1]
            if not str_id2 in dicc:
                untrans+=[str_id2]
            if str_id1 in dicc and str_id2 in dicc:
                uni_ids1=dicc[str_id1]
                uni_ids2=dicc[str_id2]
                ft.write('{}\t{}\t{}\n'\
                .format(uni_ids1,uni_ids2,int_score))
            elif force_write:
                if str_id1 in dicc:
                    str_id1=dicc[str_id1]
                if str_id2 in dicc:
                    str_id2=dicc[str_id2]
                ft.write('{}\t{}\t{}\n'\
                .format(str_id1,str_id2,int_score))
                
            
                
    untrans=set(untrans)
    for each in untrans:
        print(each, 'not found')
     #   fu.write(each+'\n')
    ft.close()
    #fu.close()
    

if __name__=="__main__":
    interaction_file=argv[1]
    
    file_dic=argv[2]
    if len(argv)>3:
        extra=argv[3]
    else:
        extra=None
        
    dicc=build_dict(file_dic)
    #print(len(dicc))
    write_interactions(interaction_file,dicc)
    
