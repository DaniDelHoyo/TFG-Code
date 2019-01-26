#!/usr/bin/env python3

'''
Author: Daniel Del Hoyo Gomez
Script that for each predicted effector with several targets,
create x random distributions of the same number of targets and
calculates a p-value for the number of times that these targets are
closer in the interactome. 
Different modes for treatment of targets in 
different subgraphs.

Input: python2 proximity_targets.py <Host_interactome> <Predicted_PPIs>
                <simul_num>
'''

from sys import argv
import subprocess as sbp
import networkx as nx
import random
import scipy.stats

def parse_effectors(ppis_file,eff_col=1):
    '''Return a dictionary {effector:[targets]}
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
    
def parse_all_prots(inter_file):
    '''Return list of all the proteins contained in a interaction file
    '''
    listi=[]
    with open(inter_file) as filex:
        for line in filex:
            listi+=[line.split()[0],line.split()[1]]
    return listi

def remove_zeros(all_pv,simu_times):
    '''Return a list of pvalues substituting zeros by very low values
    '''
    low_v=0.1/simu_times
    for i in range(len(all_pv)):
        if all_pv[i]==0:
            all_pv[i]=low_v
    return all_pv

def measure_distance(targets,inter_file):
    '''Return list with distance in interactome of all pairs of targets
    
    targets: list with target proteins
    inter_file: string, file containing the host interactome
    '''
    distances=[]
    graph=nx.read_weighted_edgelist(inter_file)
    for i in range(len(targets)-1):
        for j in range(i+1,len(targets)):
            try:
                distances+=[nx.shortest_path_length(graph,targets[i],targets[j])]
            except:
                distances+=['No']
    return distances

def set_numeric_nopath(distances,value):
    '''Set an arbitrary value of distance when it is 'no' (no path)
    '''
    for n,i in enumerate(distances):
        if i=='No':
            distances[n]=value
    return distances

def maximum_distance(inter_file):
    '''Return a int with the maximum distance between two nodes in all
    the graph'''
    graph=nx.read_weighted_edgelist(inter_file)
    eccen={}
    for g in nx.connected_component_subgraphs(graph):
        eccen.update(nx.eccentricity(g))
    
    mx=max(eccen.values())
    return mx

def comp_distances(dist_real,dist_simu,nop_value=10):
    '''Return a boolean whether the real distance is smaller than the simulated
    '''
    dist_real=set_numeric_nopath(dist_real,nop_value)
    dist_simu=set_numeric_nopath(dist_simu,nop_value)
    
    if sum(dist_real)<sum(dist_simu):
        return True
    else:
        return False

def write_pvals(p_vals,outname):
    with open(outname,'w') as f:
        listi=list(p_vals)
        listi.sort()
        for p in listi:
            f.write('{}\t{}\n'.format(p,p_vals[p]))

if __name__=="__main__":
    inter_file=argv[1]
    pred_ppis=argv[2]
    tot_reps=int(argv[3])
    outname='proximity_'+pred_ppis.split('/')[-1]
    
    subgraph_dist=maximum_distance(inter_file)+1
    print('Simulation of {} random dispositions of targets for each effector'.format(tot_reps))
    print('Distance for targets in different subgraphs: '+str(subgraph_dist))
    
    all_host_prots=parse_all_prots(inter_file)
    eff_tar_dic=parse_effectors(pred_ppis)
    
    random.seed(0)
    p_vals={}
    for effector in eff_tar_dic:
        amount_targets=len(eff_tar_dic[effector])
        if amount_targets>1:
            targets=eff_tar_dic[effector]
            real_dist=measure_distance(targets,inter_file)
            
            #Create random distributions of targets
            better_random=0
            for rep in range(tot_reps):
                simu_tagets=random.sample(all_host_prots,amount_targets)
                simu_dist=measure_distance(simu_tagets,inter_file)
                
                if not comp_distances(real_dist,simu_dist,subgraph_dist):
                    better_random+=1
            
            p_vals[effector]=better_random/float(tot_reps)
            print('{}\t{}'.format(effector,p_vals[effector]))
            #print(' '.join(['\t']+targets+['\n']))
            
    
    
    all_pv=p_vals.values()
    all_pv=remove_zeros(all_pv,tot_reps)
    sta,final_pv=scipy.stats.combine_pvalues(all_pv)
    print(sta,final_pv)
    p_vals['Fisher Combined']=final_pv
    write_pvals(p_vals,outname)
    
