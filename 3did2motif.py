#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that converts the original 3did motif domain interactions 
in motif_name   regular expresion
Input: python3 3did2motif.py <DMI.txt> <out_name>
'''

from sys import argv

def construct(ifile,ofile):
    with open(ifile,"r") as filex:
        f=open(ofile,"w")
        for line in filex:
            if line.startswith("#=ID"):
                name=line.split()[3]
            elif line.startswith("#=PT"):
                pattern=line.split()[1]
                f.write(name+'\t'+pattern+"\n")
        f.close()
    return

if __name__=="__main__":
    ifile=argv[1]
    if len(argv)>2:
        ofile=argv[2]
    else:
        ofile="3did2motif.txt"
    
    construct(ifile,ofile)
    
