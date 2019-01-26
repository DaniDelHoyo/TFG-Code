#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that check the scores of known interactions in predicted 
interaction files (dyer and xiaopan). Compare them to the average of the
scores in the predicted interactome
Input: python3 compare_known.py <known_interactions> outname
        <pred_inter1> <pred_inter2> ....
'''

from sys import argv,exit
import subprocess as sbp

def parse_known(kn_file):
    '''Return a list [host_prot|path_prot,]
    '''
    listi=[]
    with open(kn_file) as filex:
        for line in filex:
            line=line.split()
            listi.append('|'.join(line))
    return listi

def parse_predict(pred_file):
    '''Return a dictionary of dictionaries {host_prot:{path_prot:score}}
    also with: {average:average_score}
    '''
    dic={}
    con,sctot=0,0
    with open(pred_file) as filex:
        for line in filex:
            line=line.split()
            hp,pp,sc=line
            if hp in dic:
                dic[hp][pp]=sc
            else:
                dic[hp]={pp:sc}
            con+=1
            sctot+=float(sc)
    dic['average']=str(round(sctot/con,6))
    return dic

def write_output(kn_list,cubic_dic,meth_list,outname):
    '''Write a file in form
    Prot_pair   Dyer_score  Dyer_avg    Xiaopan_score   Xiaopan_avg
    '''
    f_tmp=open(outname+'tmp.txt','w')
    with open(outname+'tmp2.txt','w') as f:
        avgs={}
        f.write('Prot_pair')
        for each in meth_list:
            f.write('\t'+each.split('/')[-1])
            avgs[each]=cubic_dic[each]['average']
        f.write('\n')
        f.write('Averages_set')
        for av in meth_list:
            f.write('\t'+avgs[av])
        f.write('\n')
        
        con=0
        avgscores=[0]*len(meth_list)
        for pair in kn_list:
            hp,pp=pair.split('|')[:2]
            scores=['Non']*len(meth_list)
            i=0
            write=False
            for meth in meth_list:
                if hp in cubic_dic[meth] and pp in cubic_dic[meth][hp]:
                    scores[i]=round(float(cubic_dic[meth][hp][pp]),4)
                    avgscores[i]+=round(float(cubic_dic[meth][hp][pp]),4)
                    write=True
                    con+=1
                i+=1
                    
            
            if write:
                f_tmp.write(pair)
                for sc in scores:
                    f_tmp.write('\t{}'.format(sc))
                f_tmp.write('\n')
        
        f_tmp.close()
        f.write('Averages_known')
        for sc in avgscores:
            f.write('\t'+str(round(sc/con,6)))
        f.write('\n')
        f.close()
        sbp.check_call('cat {} {} > {}'\
        .format(outname+'tmp2.txt',outname+'tmp.txt',outname),shell=True)
        sbp.check_call('rm {} {}'\
        .format(outname+'tmp.txt',outname+'tmp2.txt'),shell=True)




if __name__=="__main__":
    kn_file=argv[1]
    outname=argv[2]
    cubic_dic,meth_list={},[]
    for i in range(3,len(argv)):
        meth_list+=[argv[i]]
        cubic_dic[argv[i]]=parse_predict(argv[i])
    
    kn_list=list(set(parse_known(kn_file)))
    
    write_output(kn_list,cubic_dic,meth_list,outname)
    
