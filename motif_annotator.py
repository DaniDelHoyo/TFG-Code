#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that annotate a protein fasta file with motifs present in a 
file. If a NetSurfP output file is provided, motifs are filtered
to >=50% exposure
Input: python3 motif_annotator.py <fasta_file> <Motifs_file> 
        <NetSurfP_File>|0 <Ocurrence_threshold>|0
#If output already exists, able to filter that output file
Input(2): python3 motif_annotator.py <motif_annotation_output> 
        <Motifs_file> <NetSurfP_File>|0 <Ocurrence_threshold>|0
NetSurfP and ocurrence filters could be skipped using 0 as argument
'''

import re
from sys import argv,exit
import hashlib
import subprocess as sbp

def parse_net(net_file):
    '''Parse the netsurfp output to make a dictionary with {prot:EBEBE}
    Being B=Buried and E=Exposed
    
    net_file: string, name of output file of netsurfp
    '''
    with open(net_file) as filex:
        exp_dic={}
        for line in filex:
            if line.startswith('#'):
                pass
            else:
                line=line.split()
                #Define protein name
                if line[0]!='X':
                    prot_name=line[2].strip()
                    
                else:
                    prot_name=line[1].strip()
                #Interproscan does extrange things when '_' or '|' is in
                #the descriptor, this is a fixing for the ap_known set
                #if '_' in prot_name:
                 #   prot_name=prot_name.split('_')[1]
                #Define exposed or buried (if unknown X)
                if prot_name in exp_dic:
                    exp_dic[prot_name]+=line[0]
                else:
                    exp_dic[prot_name]=line[0]
    return exp_dic

def make_motif_dic(mot_file):
    with open(mot_file) as filex:
        filex.readline()
        mot_dic={}
        for line in filex:
            line=line.split()
            mot_dic[line[0]]=line[1]
    return mot_dic

def parse_fasta(prot_file):
    '''Returns a dict of strings: sequences

    prot_file: string, filename of fasta'''
    seqs={}
    with open(prot_file) as filex:
        for line in filex:
            if line.startswith('>'):
                descr=line[1:].strip()
                seqs[descr]=''
            else:
                seqs[descr]+=line.strip().upper()
    return seqs

def annotator(matchy,seq,iseq,imotif,out_name):
    '''Create a file with annotated sequences with motifs,
    similar format than Interproscan output tsv'''
    with open(out_name,"a") as f:
        md5=hashlib.md5(seq.encode('utf-8')).hexdigest()
        #Another fixing for ap_known set (uniprot)
        #if '|' in iseq:
        #    iseq=iseq.split('|')[1]
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'\
        .format(iseq,md5,len(seq),imotif,matchy.group(0),\
        matchy.start(),matchy.end()))
    

def rastreator(seqs_dic,mot_dic,out_name):
    '''Rastreate the sequences with the motifs,
    if there's match it calls the function to write'''
    for iseq in seqs_dic:
        seq=seqs_dic[iseq]
        for imotif in mot_dic:
            motif=re.compile(mot_dic[imotif])
        
            matchy=re.search(motif,seq)
            if matchy!=None:
                annotator(matchy,seq,iseq,imotif,out_name)
    return (out_name)

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
    return

def calc_exposure(exp_motif):
    '''Return a float with the proportion of exposure of a motif
    
    exp_motif: string of form EEBBBXEE
    X positions are calculated depending of adyacent positions
    '''
    cont=0
    for i in range(len(exp_motif)):
        #If the exposure is unknown (X), calculated by mean of adyacents
        if exp_motif[i]=='X':
            if i>0 and i<len(exp_motif)-1:
                ant=exp_motif[i-1]
                post=exp_motif[i+1]
            else:
                if i>0:
                    ant=exp_motif[i-1]
                    post=exp_motif[i-1]
                else:
                    ant=exp_motif[i+1]
                    post=exp_motif[i+1]
            #0.25 if B adyacents, 0.75 if E adyacents, 0.5 if both
            if ant!=post:
                cont+=0.5
            elif ant=='E':
                cont+=0.75
            elif ant=='B':
                cont+=0.25
                
        elif exp_motif[i]=='E':
            cont+=1
    
    return(cont/len(exp_motif))
    
def surface_filter(out_name,exp_dic):
    '''Filter the motifs annotation by asking the motif to be at least 
    50% exposed based on netsurfp analysis
    
    out_name: string, filename of motif annotation without filtering
    exp_dic: dictionary, {seq_name:EBEBEBEBX}
    '''
    out_name2=out_name.replace('.tsv','_filt.tsv')
    f=open(out_name2,"w")
    with open(out_name) as filex:
        for linel in filex:
            line=linel.split()
            prot_name=line[0].strip()
            if prot_name in exp_dic:
                pass
            else:
                prot_name=prot_name.replace('|','_')
            exp_motif=exp_dic[prot_name][int(line[-2]):int(line[-1])]
            if calc_exposure(exp_motif)>=0.5:
                f.write(linel)
    return (out_name2)

def calc_max_ocur(thres,ocur_list,ocurren):
    '''Return the maximum value of ocurrences of a motif to get a 
    threshold given the quantile
    
    thres: float, quantile
    ocur_list: list, ocurrences of motifs ordered
    '''
    if ocurren:
        #Quantile on annotations (take only motifs for 0.75 first annotations)
        max_mot=thres*sum(ocur_list)
        i=0
        while max_mot>0:
            max_mot-=ocur_list[i]
            i+=1
    else:
        #Quantile on motifs (take only first 0.75 motifs)
        i=int(thres*len(ocur_list))
        
    return ocur_list[i]

def write_ocur_filt(good_mot,ifile,ofile):
    '''Filter the tsv annotation file writting only motifs given in a 
    list
    good_mot: list of strings, motifs to write
    ifile: string, file without filtering
    ofile: string, output file
    '''
    f=open(ofile,'w')
    with open(ifile) as filex:
        for line in filex:
            if line.split()[3] in good_mot:
                f.write(line)
    f.close()

def ocurrences_filter(thres,ifile,num_of_proteins):
    '''Filter the motif annotation retaining only those motifs that are
    annotated less than a threshold
    
    thres: float/int, if less than 1 is a quantile, if more is the 
        maximum number of occurrences
    ifile: string, file where the motif annotation is
    '''
    #Getting motif histogram in ordered lists
    tmpi='tmp'+ifile
    cmd='cut -f 4 {} |sort|uniq -c| sort -n > {}'.format(ifile,tmpi)
    xd=sbp.check_call(cmd,shell=True)
    ocur_list,mot_list=[],[]
    with open(tmpi) as filex:
        for line in filex:
            ocur_list+=[int(line.split()[0])]
            mot_list+=[line.split()[1]]
    sbp.check_call('rm '+tmpi,shell=True)
    
    #If threshold given as float, the maximum number of ocurrences is
    #the fraction of the total number of sequences
    if thres<1:
        thres=round(thres*num_of_proteins)
    
    print('Filtering by number of occurrences ({})'\
        .format(thres))
    #Choosing motifs only with lower number of ocurrences
    idx_until=0
    while ocur_list[idx_until]<=thres:
        idx_until+=1
    good_motifs=mot_list[:idx_until]  
    
    #Writing file only with good motifs
    ofile=ifile.replace('.tsv','_{}.tsv'.format(int(thres)))
    write_ocur_filt(good_motifs,ifile,ofile)
    
    return ofile
    
if __name__=="__main__":
    #Building and naming files
    prot_file=argv[1]
    mot_file=argv[2]
    net_file=argv[3]
    ocur_thres=argv[4]
    
    if not '.tsv' in prot_file:
        if len(argv)>5:
            out_name=argv[5]
        else:
            out_name=prot_file.split("/")[-1]
            if '.fasta' in out_name:
                out_name=out_name.replace('.fasta','_motifs.fasta.tsv')
            elif '.fa' in out_name:
                out_name=out_name.replace('.fa','_motifs.fa.tsv')
            else:
                exit('Incorrect fasta file')
                
        with open(out_name,"w") as f:
            pass
        #Building dictionaries
        mot_dic=make_motif_dic(mot_file)
        seqs_dic=parse_fasta(prot_file)
        
        #Searching and annotating motifs in sequences
        print('Running total annotation')
        out_name=rastreator(seqs_dic,mot_dic,out_name)
    else:
        out_name=prot_file
    cmd='cut -f 1 {} | sort | uniq | wc -l'.format(out_name)
    num_of_proteins=int(sbp.check_output(cmd,shell=True).decode())
    print(num_of_proteins,'proteins annotated')
    #Filtering by netsurfp exposure of residues
    if net_file!='0':
        print('Filtering by surface exposure')
        exp_dic=parse_net(net_file)
        
        out_name=surface_filter(out_name,exp_dic)
    #Filtering by number of ocurrences
    if ocur_thres!='0':
        out_name=ocurrences_filter(float(ocur_thres),out_name,num_of_proteins)
    
    #Creating and histogram with numbers of motif per sequence
    #tot_numb=total_proteins(prot_file)
    histo_dic=get_histogram(out_name,num_of_proteins)
    write_csv(histo_dic,out_name.replace('.tsv','_histo.csv'))
        
