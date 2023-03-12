import os
import time
import pandas as pd
import json
import networkx as nx
import numpy as np

import scipy
from scipy import stats

from multiprocessing.pool import ThreadPool
from collections import Counter

import sys

# LOAD MY DBs

path_to_dbuilder = #path to dbuilder

sys.path.append(path_to_dbuilder)

from src import extract_uniref    as uniref
from src import extract_uniprot   as uniprot
from src import extract_uniparc   as uniparc

target_uniref = 'UniRef50'

MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077

uniref_db = uniref.uniref_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT, name = target_uniref)
uniref_db.index_db()

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniprot_db.index_db()

uniparc_db = uniparc.uniparc_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniparc_db.index_db()


# GET INPUTS
AFDB_data = sys.argv[1]
uniprot_tax = sys.argv[2]
threads   = int(sys.argv[3])
infolder  = AFDB_data.split('/')[-2]

# LOAD INPUTS

print('1. Loading AFDB data')
AFDB90_CC = pd.read_csv(AFDB_data, dtype = {'communityID': str})
AFDB90_CC = AFDB90_CC.sort_values(by='unirefID')
AFDB90_CC = AFDB90_CC.set_index("unirefID") 

print('2. Loading taxonomy data')
taxonomy = pd.read_csv(uniprot_tax, index_col=0)

print('3. Getting outlier data')
outliers = '/scicore/home/schwede/durair0000/projects/turtle_tools/afdb-geometricus/data/outlier_results/'
outliers_data = {}

files = sorted(os.listdir(outliers))
for i, file in enumerate(files):
    if i % 10 == 0:
        print(i, len(files))
        
    with open('{}/{}'.format(outliers, file)) as f:
        for line in f:
            line = line.strip().split()
            uniprot_id = line[0].split('-')[0]
            if uniprot_id not in outliers_data:
                outliers_data[uniprot_id] = [float(line[1])]
            else:
                outliers_data[uniprot_id].append(float(line[1]))


# DEFINE ROUTINES

def chunk_list(l, n):

    chunks = np.array_split(np.array(l), n)

    chunks = [list(chunk) for chunk in chunks]
    return chunks

def get_tax_from_menzi(subgr_members, max_chunck_size=10000):
    
    superkingdoms = []
    
    for unirefID in subgr_members.index:
        
        curr_accs = uniref_db.query(unirefID)[0]['data']['UNIREF']['UniRef100']
        curr_accs = [i.split('_')[-1] for i in curr_accs]

        if len(curr_accs) > max_chunck_size:
            ratio    = len(curr_accs)/max_chunck_size
            if ratio - round(ratio) > 0:
                n_chunks = round(ratio) + 1
            else:
                n_chunks = round(ratio)                
            chuncks  = chunk_list(curr_accs, n_chunks, counts=None)
        else:
            chuncks = [curr_accs]

        for chunck in chuncks:
            
            up_docs = uniprot_db.col.find({'_id': {'$in': chunck}})
            for doc in up_docs:
                acc = doc['_id']
                try:
                    tax = doc['data']['TAXID'][2][0]
                    superkingdoms.append(tax)
                except:
                    pass

            uparc_docs = uniparc_db.col.find({'_id': {'$in': curr_accs}})
            for doc in uparc_docs:
                acc = doc['_id']
                try:
                    tax = doc['data']['TAXID'][2][0]
                    superkingdoms.append(tax)
                except:
                    pass
        
    try:
        count = Counter(superkingdoms)    
        return count.most_common(1)[0][0], count.most_common(1)[0][1]*100/len(superkingdoms)
    except:
        return np.nan, np.nan
    
def get_comunities_summary_for_communities(arguments):
    
    target_communities = arguments[0]
    AFDB90_CC    = arguments[1]
    taxonomy     = arguments[2]
    outlier_data = arguments[3]
    thread_id    = arguments[4]
    
    curr_data = AFDB90_CC.loc[AFDB90_CC.communityID.isin(target_communities)]
    curr_tax  = taxonomy.loc[taxonomy.communityIDs.isin(target_communities)]

    start = time.time()

    communities_summary = {'Community': [], 'Subgraph': [], 'Avg_darkness': [], 'SD_darkness': [], 
                           'Avg_outlier_score': [], 'SD_outlier_score':[], 'N_members': [], 'TM': [], 
                           'SP': [], 'Median_length': [], 'MAD_length': [], 'Median_darkness': [], 'MAD_darkness': [],
                           'Median_representative': [], 'Longest_representative': [], 'Median_rep_title': [],
                           'Mode_superkingdom': [], 'Freq_superkingdom': []}

    n_expected = len(target_communities)
    for i, community_class in enumerate(target_communities):

        subgr_members = curr_data.loc[curr_data.communityID == community_class]

        communities_summary['Subgraph'].append(community_class.split('[')[0])
        communities_summary['Community'].append(community_class)

        communities_summary['Avg_darkness'].append(np.mean(subgr_members.FULL_noDUF.astype(float)))
        communities_summary['SD_darkness'].append(np.std(subgr_members.FULL_noDUF.astype(float)))

        communities_summary['Median_darkness'].append(np.median(subgr_members.FULL_noDUF.astype(float)))
        communities_summary['MAD_darkness'].append(stats.median_abs_deviation(subgr_members.FULL_noDUF.astype(float)))

        outlier_scores = [outlier_data[i] for i in subgr_members.AF2_longest_best70 if i in outlier_data]
        if len(outlier_scores) > 0:    
            communities_summary['Avg_outlier_score'].append(np.mean(outlier_scores))
            communities_summary['SD_outlier_score'].append(np.std(outlier_scores))
        else:
            communities_summary['Avg_outlier_score'].append(np.nan)
            communities_summary['SD_outlier_score'].append(np.nan)

        communities_summary['N_members'].append(len(subgr_members))

        communities_summary['TM'].append((len(subgr_members)-list(subgr_members.TM).count(0))*100/len(subgr_members))
        communities_summary['SP'].append((len(subgr_members)-list(subgr_members.SP).count(0))*100/len(subgr_members))

        median = np.median(subgr_members.AF2_longest_best70_len.astype(float))
        communities_summary['Median_length'].append(median)
        communities_summary['MAD_length'].append(stats.median_abs_deviation(subgr_members.AF2_longest_best70_len.astype(float)))

        subgr_members['dist_to_median'] = abs(subgr_members.AF2_REP_best_len - median) 
        
        median_rep = subgr_members.sort_values(by='dist_to_median', ascending=True).AF2_longest_best70[0]
        communities_summary['Median_representative'].append(median_rep)
        
        longest_rep = subgr_members.sort_values(by='AF2_REP_best_len', ascending=False).AF2_longest_best70[0]
        communities_summary['Longest_representative'].append(longest_rep)
        
        median_title = uniprot_db.query(median_rep)[0]['data']['NAME']['TITLE']
        communities_summary['Median_rep_title'].append(median_title)
        
        try:
            subgrp_tax = curr_tax.loc[curr_tax.communityIDs == community_class]
            tax = subgrp_tax.superkingdom.mode()[0]
            frq = subgrp_tax['superkingdom'].value_counts()[tax]*100/len(subgrp_tax)
        except:
            tax, frq = get_tax_from_menzi(subgr_members)
            
        communities_summary['Mode_superkingdom'].append(tax)
        communities_summary['Freq_superkingdom'].append(frq)
        
        if i % 100 == 0:
            numb_seconds = time.time() - start
            time_to_end = round(((numb_seconds/(i+1))*n_expected)-numb_seconds)
            print('thread {}:'.format(thread_id), i+1, n_expected, 'CURR COMMUNITY:', community_class, 'CURR TITLE', median_title, ' ... Time passed: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))), 'Expected to finish in: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(time_to_end)))-1, int(time.strftime('%d', time.gmtime(time_to_end)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(time_to_end))), flush = True)  

    return communities_summary
    

def get_comunities_summary(AFDB90_CC, outlier_data, taxonomy, threads):
    
    target_communities = list(set(AFDB90_CC.communityID))
    separated_jobs = chunk_list(target_communities, threads)
    
    list_arguments = [i for i in zip(separated_jobs, [AFDB90_CC for job in separated_jobs], [taxonomy for job in separated_jobs], [outlier_data for job in separated_jobs], range(threads))]
    
    pool = ThreadPool(threads)
    results = pool.imap_unordered(get_comunities_summary_for_communities, list_arguments)
    
    all_results = {}
    for dic in results:
        for key in dic:
            if key not in all_results:
                all_results[key] = dic[key]
            else:
                all_results[key] += dic[key]
    
    all_results = pd.DataFrame(all_results)
    all_results = all_results.set_index('Community')
    all_results = all_results.sort_values(by='N_members', ascending=False)
    
    return all_results  
    
    
# GET COMMUNITIES SUMMARY

print('3. Getting communities summary')

if not os.path.isfile('{}/communities_summary_noreps.csv'.format(infolder)):
    communities_summary = get_comunities_summary(AFDB90_CC, outliers_data, taxonomy, threads=threads)
    communities_summary.to_csv('{}/communities_summary_noreps.csv'.format(infolder))
    
else:
    communities_summary = pd.read_csv('{}/communities_summary_noreps.csv'.format(infolder), dtype = {'Community': str})
    communities_summary = communities_summary.set_index("Community")  

communities_summary.to_csv('{}/communities_summary.csv'.format(infolder))