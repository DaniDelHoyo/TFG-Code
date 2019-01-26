#!/usr/bin/env python3

from sys import argv

'''
Author: Daniel Del Hoyo Gomez
Script that analyzes the frequency of a motif annotation
Input: python3 motif_frequency.py <motif_anotated tsv> 
'''

def motif_frequency(filename):
    '''Return a dictionary with the frequency of appearence of a motif
    in a set of sequences. 
    
    filename: string, file resulting of motif annotation
    '''
    freq_dic={}
    len_dic={}
    with open(filename) as filex:
        for line in filex:
            line=line.split()
            if line[3] in freq_dic:
                freq_dic[line[3]]+=1
            else:
                freq_dic[line[3]]=1
                
            if len(line[4]) in len_dic:
                len_dic[len(line[4])]+=1
            else:
                len_dic[len(line[4])]=1
                
    return freq_dic,len_dic

def write_freq_file(freq_dic,len_dic,out_name='motif_freq.csv'):
    '''Write a csv file with the frequency of each motif ina  set
    Also writes the histogram of frequencies
    
    freq_dic: dictionary, contains {motif:appereances}
    out_name: string, output filename
    '''
    maxi_freq=max(list(freq_dic.values()))
    hist_dic=[0]*maxi_freq
    maxi_len=max(list(len_dic.keys()))
    hist_len=[0]*maxi_len
    for lens in len_dic:
        hist_len[lens-1]=len_dic[lens]
        
    towrite=[]
    for motif in freq_dic:
        towrite+=['{},{},,'.format(motif,freq_dic[motif])]
        hist_dic[freq_dic[motif]-1]+=1
    
    f=open(out_name,'w')
    for i in range(max(maxi_freq,len(towrite))):
        if i<len(towrite):
            f.write(towrite[i])
        else:
            f.write(',,,')
            
        if i<len(hist_dic):
            f.write('{},{},,'.format(i+1,hist_dic[i]))
        else:
            f.write(',,,')
            
        if i<len(hist_len):
            f.write('{},{}\n'.format(i+1,hist_len[i]))
        else:
            f.write('\n')
            
    f.close()
    
    
if __name__=="__main__":
    ifile=argv[1]
    freq_dic,len_dic=motif_frequency(ifile)
    outname=ifile.replace('.tsv','')+'_freq.csv'
    write_freq_file(freq_dic,len_dic,outname)
