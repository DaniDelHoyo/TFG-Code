#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Student number: 960824-177-040
Script that annotates a interactome file with a gene_association file 
from GOA uniprot
Input is:
    python3 go_enrichment.py <uniprot_host_interactome> <gaf_file> 
            <uniprot_pred_interactome> <interesting_gos>
'''

from sys import argv
import subprocess as sbp

def parse_infect_gos(infection_gos):
    '''Return a list with gos in first column'''
    listi=[]
    with open(infection_gos) as filex:
        for line in filex:
            listi+=[line.split()[0]]
    return listi

def parse_gobranches(b_file):
    '''Return a dictionary {GO:branch}
    '''
    dic={}
    with open(b_file) as filex:
        for line in filex:
            line=line.split()
            dic[line[0]]=line[1]
    return dic
    
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

def parse_inter(inter_file,col=None):
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
    '''Returns a go dictionary only with proteins in list inter_prots
    {prot:[GOs]}
    '''
    dic={}
    for prot in inter_prots:
        if prot in go_dic:
            dic[prot]=go_dic[prot]
    return dic

def calc_gos_frequency(prot_dic):
    '''Return a dictionary with {GO1:proportion of proteins with GO1}
    '''
    total_prots=len(prot_dic)
    gof_dic={}
    for prot in prot_dic:
        for go in prot_dic[prot]:
            if go in gof_dic:
                gof_dic[go]+=1
            else:
                gof_dic[go]=1
    for go in gof_dic:
        gof_dic[go]=gof_dic[go]/total_prots
    return gof_dic

def write_dic(dic,outname,good_gos):
    '''Write a tsv file with key value of a dic
    '''
    with open(outname,'w') as f:
        for go in dic:
            if go in good_gos:
                tiki='Inf_rel'
            else:
                tiki='Non_rel'
            f.write('{}\t{}\t{}\t{}\n'\
            .format(go,dic[go][0],str(dic[go][1]),tiki))
    sbp.check_call('sort -k 3 -g {} > pr'.format(outname),shell=True)
    sbp.check_call('mv pr {}'.format(outname),shell=True)        

def write_stats(tp,fp,good_gos,gos_present,outname):
    '''Write stats on a file
    '''
    good_present=set(good_gos)&set(gos_present)
    with open(outname,'w') as f:
        f.write('GOs enriched\n')
        f.write('Infection related GOs: {}/{}/{}\n'\
        .format(tp,len(good_present),len(good_gos)))
        noninf_len=len(gos_present)-len(good_present)
        f.write('Non infection related GOs: {}/{}/{}\n'\
        .format(fp,noninf_len,len(gos_present)))
        f.write('\nPrecision: {}\n'.format(round(tp/(tp+fp),3)))
        f.write('Recall: {}\n'.format(round(tp/(len(good_present)),3)))
        f.write('Accuracy: {}\n'.format(round((tp+noninf_len-fp)/len(gos_present),3)))


if __name__=="__main__":
    inter=argv[1]
    gaf=argv[2]
    targets_file=argv[3]
    infect_gos=argv[4]
        
    branches=[]
    i=0
    while len(infect_gos.split('_')[i])==2:
        branches+=[infect_gos.split('_')[i]]
        i+=1
        
    if len(argv)>5:
        outname=argv[5]
    else:
        outname='goenrich_'+targets_file.split('/')[-1]
        if 'chil' in infect_gos:
            outname='chil_'+outname
        outname='_'.join(branches+[outname])
        print(outname)
    
    for i in range(len(branches)):
        if branches[i]=='bp':
            branches[i]='biological_process'
        elif branches[i]=='mf':
            branches[i]='molecular_function'
        elif branches[i]=='cc':
            branches[i]='celular_component'
    
    #Parsing GO term with GO function and branches
    go_func=parse_gofun('GO_functions')
    go_branches=parse_gobranches('go_links/GO_branches.tsv')
    good_gos=parse_infect_gos(infect_gos)
    
    #Taking interactome and targeted proteins lists
    inter_prots=parse_inter(inter)
    target_prots=parse_inter(targets_file,col=0)
    
    #Taking a dictionary with {protein:[gos]}
    go_dic=parse_gaf(gaf)
    #Taking a GO dictionary only with proteins in interactome,targets
    inter_prot_dic=filter_dic(go_dic,inter_prots)
    target_prot_dic=filter_dic(go_dic,target_prots)
    
    #List of GO terms in targets / interactome
    #Enrichment in targeted subset of infection related GOs
    gos_targeted,hitS=[],0
    for prot in target_prot_dic:
        for go in target_prot_dic[prot]:
            gos_targeted+=[go]
            if go in good_gos:
                hitS+=1
    
    gos_interactome,hitP=[],0
    for prot in inter_prot_dic:
        for go in inter_prot_dic[prot]:
            gos_interactome+=[go]
            if go in good_gos:
                hitP+=1
    gos_present=list(set(gos_interactome))
    
    totalSize=len(gos_interactome)
    sampleSize=len(gos_targeted)
    cmd='Rscript go_enrich.R {} {} {} {}'\
    .format(hitS,hitP,totalSize,sampleSize)
    p_val=sbp.check_output(cmd,shell=True).decode().split()[1]
    print('Infection related GOs in targets: {}/{}'\
    .format(hitS,sampleSize))
    print('Infection related GOs in interactome: {}/{}'\
    .format(hitP,totalSize))
    print('Fisher pvalue enrichment: {}'.format(p_val))
    
    #Enrichment for each GO term present in targets
    gos_enriched={}
    tp,fp=0,0
    for go in set(gos_targeted):
        if go_branches[go] in branches:
            hitS=gos_targeted.count(go)
            hitP=gos_interactome.count(go)
            cmd='Rscript go_enrich.R {} {} {} {}'\
            .format(hitS,hitP,totalSize,sampleSize)
            p_val=sbp.check_output(cmd,shell=True).decode().split()[1]
            if float(p_val)<0.1:
                try:
                    fun=go_func[go]
                except:
                    fun='???'
                gos_enriched[go]=[fun,p_val]
                if go in good_gos:
                    tp+=1
                else:
                    fp+=1
                    
    write_dic(gos_enriched,outname,good_gos)
    
    good_present=set(good_gos)&set(gos_present)
    noninf_len=len(gos_present)-len(good_present)
    tn=noninf_len-fp
    print('\nGOs enriched')
    print('Infection related GOs: {}/{}/{}'\
    .format(tp,len(good_present),len(good_gos)))
    
    print('Non infection related GOs: {}/{}/{}'\
    .format(fp,noninf_len,len(gos_present)))
    print('\nPrecision: {}'.format(round(tp/(tp+fp),3)))
    print('Recall: {}'.format(round(tp/(len(good_present)),3)))
    print('Specificity: {}'.format(round(tn/(fp+tn),3)))
    print('Accuracy: {}'.format(round((tp+tn)/len(gos_present),3)))
    
    write_stats(tp,fp,good_gos,gos_present,outname.replace('.tsv','_stats.tsv'))




