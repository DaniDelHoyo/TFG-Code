#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that 

Input: python3 correlation_common_effectors.py <predicted_interactions>
        <pathogen_expression_file> <host_expression_file> 
        <simulation_times>
'''

from sys import argv
import subprocess as sbp
import scipy.stats
import random
from numpy import isnan

def parse_expression(exp_file):
    '''Return a dictionary {gene_id:[expression_counts]}
    '''
    dic={}
    with open(exp_file) as filex:
        for line in filex:
            if not line.startswith('gene'):
                line=line.split()
                dic[line[0]]=line[1:]
    return dic

def parse_effectors(ppis_file,eff_col=0):
    '''Return a dictionary {target:[effectors]}
    '''
    dic={}
    #If 0->1. If 1->0
    tar_col=abs(eff_col-1)
    with open(ppis_file) as filex:
        for line in filex:
            line=line.split()
            if line[eff_col] in dic:
                dic[line[eff_col]]+=[line[tar_col]]
            else:
                dic[line[eff_col]]=[line[tar_col]]
    return dic

def random_correlation_pvalue(effs,pval,path_exp_dic,times=100):
    '''Pick a random effector and perform the Spearman correlation with 
    each of the known effectors.
    Compare the p-value of the random Spearman with the real case.
    Return a p-value for the probability of finding a better correlation
    randomly
    '''
    random.seed(0)
    r_effectors=random.sample(list(path_exp_dic.keys()),times+1)
    results=[]
    for effi in effs:
        expi=path_exp_dic[effi]
        #The analysis will be done times, the sampling is times+1 because 
        #the effi itself may be sampled
        resul,con=0,0
        while con<=times:
            if r_effectors[con]==effi:
                con+=1
            effx=r_effectors[con]
            #print(effi,effx)
            expx=path_exp_dic[effx]
            r_pval=scipy.stats.spearmanr(expi,expx)[1]
            #print(pval,r_pval)
            
            if r_pval<=pval:
                resul+=1
            con+=1
        results+=[resul/times]
    return sum(results)-results[0]*results[1]

def rho_distribution(host_exp_dic,path_exp_dic,times=1000):
    '''Return a list with absolute values of Spearman rhos for random
    pairs of genes
    
    path_exp_dic: dictionary {gene:[expressions]}
    '''
    #Made a list and sorted for consistency with seed
    path_prots=list(path_exp_dic)
    path_prots.sort()
    host_prots=list(host_exp_dic)
    host_prots.sort()
    random.seed(0)
        
    rho_dis=[]
    for t in range(times):
        #Chosing a random pair
        while True:
            i,j=random.sample(host_prots,1)[0],random.sample(path_prots,1)[0]
            if i!=j:
                expi,expj=host_exp_dic[i],path_exp_dic[j]
                expi,expj=list(map(float,expi)),list(map(float,expj))
                rho,pval=scipy.stats.pearsonr(expi,expj)
                if not isnan(rho):
                    break
            
        rho_dis+=[abs(rho)]
        
    rho_dis.sort()
    return rho_dis

def adjs_pval(rho_dis,rho):
    '''Return a float with the proportion of rho values bigger than rho
    '''
    rho=abs(rho)
    i=0
    while i<len(rho_dis) and rho>rho_dis[i]:
        i+=1
        
    return (len(rho_dis)-i)/len(rho_dis)

def remove_zeros(all_pv,simu_times):
    '''Return a list of pvalues substituting zeros by very low values
    '''
    low_v=0.1/simu_times
    for i in range(len(all_pv)):
        if all_pv[i]==0:
            all_pv[i]=low_v
    return all_pv

def write_results(spear,outname):
    '''Write a file with the results
    '''
    with open(outname,'w') as f:
        f.write('Protein_1,Protein_2,Pearson_coef,Pearson_Pvalue,Adjusted_Pvalue\n')
        listi=list(spear)
        listi.sort()
        for target in listi:
            
            for pair in spear[target]:
                pair=tuple(pair[:2])+tuple(map(lambda n: "%.4g" % n,pair[2:]))
                f.write('{},{},{},{},{}\n'\
                .format(pair[0],pair[1],pair[2],pair[3],pair[4]))
            f.write('\n')

if __name__=="__main__":
    pred_ppis_file=argv[1]
    path_exp_file=argv[2]
    host_exp_file=argv[3]
    simu_times=int(argv[4])
    if len(argv)>5:
        outname=argv[5]
    else:
        outname='expression_'+pred_ppis_file.split('/')[-1]
        outname=outname.replace('.tsv','.csv')
    
    targ_dic=parse_effectors(pred_ppis_file)
    path_exp_dic=parse_expression(path_exp_file)
    host_exp_dic=parse_expression(host_exp_file)

    rho_dis=rho_distribution(path_exp_dic,path_exp_dic,simu_times)
    rho_dis_te=rho_distribution(host_exp_dic,path_exp_dic,simu_times)
       
    all_pv1,all_pv2=[],[]
    spear={}    
    for targ in targ_dic:
        
        #Changing name to be able to parse expression
        t=targ.split('.')
        t_exp='.'.join(t[1:len(t)-1])
        
        spear[targ]=[]
        #Calculating correlation target-effector
        for eff in targ_dic[targ]:
            path_exp=path_exp_dic[eff]
            host_exp=host_exp_dic[t_exp]
            host_exp,path_exp=list(map(float,host_exp)),list(map(float,path_exp))
            #rho,pval=scipy.stats.spearmanr(host_exp,path_exp)
            rho,pval=scipy.stats.pearsonr(host_exp,path_exp)
            ad_pval=adjs_pval(rho_dis_te,rho)
            spear[targ]+=[(targ,eff,rho,pval,ad_pval)]
            all_pv1+=[ad_pval]
        if len(targ_dic[targ])>1:
            #Calculating correlation between effectors
            for i in range(len(targ_dic[targ])-1):
                for j in range(i+1,len(targ_dic[targ])):
                    effi=targ_dic[targ][i]
                    effj=targ_dic[targ][j]
                    
                    if effi in path_exp_dic:
                        expi=list(map(float,path_exp_dic[effi]))
                    else:
                        print(effi+' not in expression data')
                    if effj in path_exp_dic:
                        expj=list(map(float,path_exp_dic[effj]))
                    else:
                        print(effj+' not in expression data')
                    expi,expj=list(map(float,expi)),list(map(float,expj))
                    #rho,pval=scipy.stats.spearmanr(expi,expj)
                    rho,pval_pea=scipy.stats.pearsonr(expi,expj)
                    ad_pv=adjs_pval(rho_dis,rho)
                    
                    spear[targ]+=[(effi,effj,rho,pval_pea,ad_pv)]
            
                    all_pv2+=[ad_pv]
                    
                    #random_correlation_pvalue((effi,effj),pval,path_exp_dic,\
                    #simu_times))]
                    
                    #print(spear[targ])
    
    write_results(spear,outname)
    
    all_pv1=remove_zeros(all_pv1,simu_times)
    all_pv2=remove_zeros(all_pv2,simu_times)
    sta1,final_pv1=scipy.stats.combine_pvalues(all_pv1)
    print('Target-Effector',sta1,final_pv1)
    sta2,final_pv2=scipy.stats.combine_pvalues(all_pv2)
    print('Effector-Effector',sta2,final_pv2)
    
    with open(outname,'a') as f:
        f.write('Target-Effector\t'+str(final_pv1)+'\n')
        f.write('Effector-Effector\t'+str(final_pv2)+'\n')




