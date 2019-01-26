#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that creates a file with number of proteins with a number of 
domains
Input: python3 domain_annotator.py <sequences_file> <interpr_analysis[Pfam]>
'''

from sys import argv
import subprocess as sbp
import os

def run_InterProScan(sequences_file, analysis):
    '''Run the software InterProScan to annotate protein domains
    
    seq_file: string, fasta filename
    '''
    outname=sequences_file[:].split('/')[-1]
    outname+='_domains.tsv'
    
    
    if analysis[0]=='all':
        apl=''
        outname=outname.replace('domains','alldomains')
    else:
        apl=' -appl '+','.join(analysis)
    ix=2
    while os.path.exists(outname):
        #print('Not running InterProScan, file with output name exists')
        outname=outname.replace('.tsv',str(ix)+'.tsv')
        ix+=1
        
    cmd='./interproscan.sh -i {}{} -f tsv -o {}'\
    .format(sequences_file,apl,outname)
    input(outname)
    sbp.check_call(cmd,shell=True)
    return outname
    

def total_proteins(sequences_file):
    '''Return the number of proteins in a fasta file
    
    sequences_file: string, filename of fasta file
    '''
    num=int(sbp.check_output('grep ">" {} | wc -l'\
    .format(sequences_file),shell=True))
    
    return num
    
def get_histogram(interpro_file,total_prots):
    '''Return a dictionary with {number_of_domains:number_of_proteins}
    
    interpro_file: string, filename of interpro output
    total_prots: int, total number of proteins in set
    '''
    dicc={}
    nozero=sbp.check_output('cat {} | cut -f1 | uniq -c | wc -l'\
    .format(interpro_file),shell=True)
    dicc[0]=total_prots-int(nozero)
    
    histo=sbp.check_output('cat {} | cut -f1 | uniq -c | sort \
    > histo.txt'.format(interpro_file),shell=True)
    
    with open('histo.txt') as filex:
        for line in filex:
            doms=int(line.split()[0])
            if doms in dicc:
                dicc[doms]+=1
            else:
                dicc[doms]=1
    
    sbp.check_call('rm histo.txt',shell=True)
    
    return dicc

def write_csv(dicc,output):
    '''Write a csv file from a dictionary: key,value
    
    dicc: dictionary
    '''
    doms=0
    with open (output,'w') as f:
        nums=list(dicc.keys())
        nums.sort()
        while doms<=nums[-1]:
            if not doms in dicc:
                f.write('{},0\n'.format(str(doms)))
            else:
                f.write('{},{}\n'.format(str(doms),str(dicc[doms])))
            doms+=1


if __name__=="__main__":
    seq_file=argv[1]
    if len(argv)>2:
        analysis=argv[2:]
    else:
        analysis=['Pfam']
    #analysis=['all']
    interpro_file=run_InterProScan(seq_file,analysis)
    out_name=interpro_file.replace('.tsv','_histo.csv')
        
    total_prots=total_proteins(seq_file)
    dicti=get_histogram(interpro_file,total_prots)
    write_csv(dicti,out_name)
        
        
        
        
        
