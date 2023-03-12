import os
import time
import pandas as pd
import json
import itertools
import networkx as nx
import numpy as np

import scipy
from scipy import stats

from ete3 import NCBITaxa

from multiprocessing.pool import ThreadPool

import sys

import warnings
warnings.filterwarnings("ignore")

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
threads   = int(sys.argv[2])
infolder  = AFDB_data.split('/')[-2]

# LOAD INPUTS

print('1. Loading AFDB data')
AFDB90_CC = pd.read_csv(AFDB_data, dtype = {'communityID': str})
AFDB90_CC = AFDB90_CC.sort_values(by='unirefID')
AFDB90_CC = AFDB90_CC.set_index("unirefID")  

# ROUTINES

def chunk_list(l, n, counts=None):
    
    print(len(l))

    print('Making chuncks')
    
    if counts is not None:
        b = [[l[i]]*counts[i] for i in range(len(l))]
        b = list(itertools.chain.from_iterable(b))
    else:
        b = l
    
    chunks = np.array_split(np.array(b), n)
#     chunks = [list(chunk) for chunk in chunks]
    
    final_chunks = []
    for i, chunk in enumerate(chunks):
        chunk = set(chunk)
        
        if i > 0:
            last_chunk = set(final_chunks[-1])
            if len(chunk.intersection(last_chunk)) > 0:
                chunk = chunk - last_chunk
            
        chunk = list(chunk)
        final_chunks.append(chunk) 
        print(len(chunk), chunk[0])
    
    sumlen = sum([len(i) for i in final_chunks])
    print(' ... Made {} chuncks ({} jobs in total)'.format(len(final_chunks), sumlen), len(l))
    
    return final_chunks

def get_taxonomy_for_unirefs(arguments,max_chunck_size = 10000):
    
    NCBI = NCBITaxa()
    
    target_unirefs = arguments[0]
    AFDB90_CC = arguments[1]
    outfolder = arguments[2]
    thread_id = arguments[3]
    
    curr_data = AFDB90_CC.loc[AFDB90_CC.index.isin(target_unirefs)]

    target_ranks = ['superkingdom','phylum','class','order','genus','species']
    
    out_json = '{}/{}_taxonomy.json'.format(outfolder, thread_id)
    out_summary = '{}/{}_taxonomy_summary.json'.format(outfolder, thread_id)
    
    try:
        taxonomy = json.load(open(out_json, 'r'))
        target_unirefs = list(set(target_unirefs)-set(taxonomy['UniRef50IDs']))
    except:
        taxonomy = {rank: [] for rank in target_ranks}
        taxonomy['uniprotIDs'] = []
        taxonomy['UniRef50IDs'] = []
        taxonomy['communityIDs'] = []

    count = 0
    n_expected = len(target_unirefs)

    start = time.time()

    for unirefID in target_unirefs:

        row = curr_data.loc[unirefID]
        curr_community = row.communityID
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
                taxid = doc['data']['TAXID'][0]

                curr_tax = {rank: np.nan for rank in target_ranks}
                try:
                    lineage = NCBI.get_lineage(taxid)
                    translation = NCBI.get_taxid_translator(lineage)
                    ranks = NCBI.get_rank(lineage)
            
                    for level in lineage:
                        if ranks[level] in target_ranks:
                            curr_tax[ranks[level]] = translation[level]
                except:
                    pass            

                for rank in curr_tax:
                    taxonomy[rank].append(curr_tax[rank])
                    
                taxonomy['uniprotIDs'].append(acc)
                taxonomy['UniRef50IDs'].append(unirefID)
                taxonomy['communityIDs'].append(curr_community)

            uparc_docs = uniparc_db.col.find({'_id': {'$in': curr_accs}})
            for doc in uparc_docs:
                acc = doc['_id']
                taxid = doc['data']['TAXID'][0]

                curr_tax = {rank: np.nan for rank in target_ranks}
                try:
                    lineage = NCBI.get_lineage(taxid)
                    translation = NCBI.get_taxid_translator(lineage)
                    ranks = NCBI.get_rank(lineage)
            
                    for level in lineage:
                        if ranks[level] in target_ranks:
                            curr_tax[ranks[level]] = translation[level]
                except:
                    pass            

                for rank in curr_tax:
                    taxonomy[rank].append(curr_tax[rank])

                taxonomy['uniprotIDs'].append(acc)
                taxonomy['UniRef50IDs'].append(unirefID)
                taxonomy['communityIDs'].append(curr_community)

            if count % 100 == 0:
                numb_seconds = time.time() - start
                time_to_end = round(((numb_seconds/(count+1))*n_expected)-numb_seconds)
                print('thread {}:'.format(thread_id), count+1, n_expected, ' ... Time passed: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))), 'Expected to finish in: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(time_to_end)))-1, int(time.strftime('%d', time.gmtime(time_to_end)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(time_to_end))), flush = True)  

            if count % 1000 == 0:
                json.dump(taxonomy, open(out_json, 'w'), indent=4)

            count+=1
    
    json.dump(taxonomy, open(out_json, 'w'), indent=4)
    
    return taxonomy

# GET TAXONOMY

separated_jobs = chunk_list(list(AFDB90_CC.index), threads, counts=list(AFDB90_CC.nUniRef100))

list_arguments = [i for i in zip(separated_jobs, [AFDB90_CC for job in separated_jobs], [infolder for job in separated_jobs], range(threads))]

pool = ThreadPool(threads)
results = pool.imap_unordered(get_taxonomy_for_unirefs, list_arguments)

all_results = {}
for dic in results:
    for key in dic:
        if key not in all_results:
            all_results[key] = dic[key]
        else:
            all_results[key] += dic[key]

taxonomy = pd.DataFrame(all_results)
taxonomy = taxonomy.set_index('uniprotIDs')

taxonomy.to_csv('{}/uniprot_community_taxonomy_map.csv'.format(infolder))

