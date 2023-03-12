path_to_dbuilder = #path to dbuilder

import sys
import os

sys.path.append(path_to_dbuilder)
from src import extract_uniref    as uniref
from src import extract_alphafold as alphafold
from src import extract_uniprot   as uniprot

import time 
import numpy as np

def write_data_to_file(data, target_uniref, outfile):
    
    if not os.path.isfile(outfile):
        with open(outfile, 'w') as outp:
            outp.write(','.join(sorted(list(data.keys()))))
            outp.write('\n')
    
    with open(outfile, 'a+') as outp:
        for i in range(len(data['unirefID'])):
            line_data = [str(data[key][i]) if data[key][i] is not None else str(np.nan) for key in sorted(list(data.keys()))]
            outp.write(','.join(line_data))
            outp.write('\n')
    
    return {key: [] for key in data} 


# LOAD TARGET UNIREF DATABASE

target_uniref = sys.argv[1] # either UniRef90 or UniRef50

MONGO_HOST = "10.1.0.202"
MONGO_PORT = 30077

uniref_db = uniref.uniref_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT, name = target_uniref)
uniref_db.index_db()

alphafold_db = alphafold.alphafold_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
alphafold_db.index_db()

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniprot_db.index_db()

# COUNT ENTRIES IN THE UNIREF DATABASE

print('\nCOUNTING DATABASE SIZE')

n_entries = uniref_db.col.count_documents({})
print(' ... FOUND {} ENTRIES'.format(n_entries))

# DEFINE OUTPUT FILE AND CHECK THE UNIREF IDS ALREADY THERE

print('\nDEFINING OUTPUT FILE')

outfile = 'data_generated/AFDBv4_pLDDT_diggestion_{}.csv'.format(target_uniref)

# GO THROUGH EACH CHUNCK AND COLLECT THE TARGET DATA

start = time.time()

print('\nGOING THROUGH THE CHUNCKS')

write_step = 50000

data = {'unirefID': [], 'median_pLDDT': [], 'max_pLDDT': [], 'min_pLDDT': [], 'delta_pLDDT': [], 'nACCs': [], 
        'nUniRef100': [], 'nUniRef90': [], 'nAF2': [],'AF2_REP_best_len': [], 'AF2_REP_worst_len': [], 
        'AF2_REP_best': [], 'AF2_REP_worst': [], 'AF2_longest_best70': [],'AF2_longest_best70_pLDDT': [], 
        'AF2_longest_best70_len': [], 'median_Evidence': []}
#         , 'uniref_rep_len': []}

curr_count = 0
for document in uniref_db.col.find():
    uniref_id = document['_id']
    uniref_dt = document['data']

    curr_count += 1

    for key in uniref_dt['DARKNESS']:                
        if key != 'pLDDTs':
            if 'AF2_REP' in key:
                if uniref_dt['DARKNESS'][key] is not None:
                    data['{}_len'.format(key)].append(uniref_dt['DARKNESS'][key]['LEN']) 
                    data[key].append(uniref_dt['DARKNESS'][key]['ACC']) 

                else:
                    data['{}_len'.format(key)].append(np.nan) 
                    data[key].append(None) 

            else:
                if key not in data:
                    data[key] = []
                try:
                    data[key].append(round(uniref_dt['DARKNESS'][key], 2))
                except:
                    data[key].append(uniref_dt['DARKNESS'][key])

        if key == 'pLDDTs':
            if len(uniref_dt['DARKNESS'][key]) > 0:
                median = np.median(uniref_dt['DARKNESS'][key])
                maximum = max(uniref_dt['DARKNESS'][key])
                minimum = min(uniref_dt['DARKNESS'][key])
                delta = min(uniref_dt['DARKNESS'][key])-max(uniref_dt['DARKNESS'][key])

                data['median_pLDDT'].append(median)
                data['max_pLDDT'].append(maximum)
                data['min_pLDDT'].append(minimum)
                data['delta_pLDDT'].append(delta)

            else:
                data['median_pLDDT'].append(np.nan)
                data['max_pLDDT'].append(np.nan)
                data['min_pLDDT'].append(np.nan)
                data['delta_pLDDT'].append(np.nan)

    # get the longest protein with a pLDDT > 70
    longest_pLDDT_best70 = None
    best_pLDDT = None
    best_n_res = 0
    
    af_docs = alphafold_db.col.find({'_id': {'$in': uniref_dt['ACC']}})
    for af_document in af_docs:
        curr_dt = af_document['data']
        
        # get plddt
        avgPLDDT = []
        n_res = 0
        for fragment in curr_dt:
            avgPLDDT.append(curr_dt[fragment]['pLDDT']['avg_pLDDT']*curr_dt[fragment]['pLDDT']['Lenght'])
            n_res += curr_dt[fragment]['pLDDT']['Lenght']

        fullprotein_pLDDT = sum(avgPLDDT)/n_res     
        
        if fullprotein_pLDDT > 70 and n_res > best_n_res:
            longest_pLDDT_best70 = af_document['_id']
            best_pLDDT = fullprotein_pLDDT
            best_n_res = n_res
    
    if best_n_res == 0:
        best_n_res = None
    
#     rep_length = len(uniprot_db.query(uniref_id.split('_')[-1])[0]['data']['SEQ'])

    # get the median evidence level
    evidence_level = []
    uniprot_docs = uniprot_db.col.find({'_id': {'$in': uniref_dt['ACC']}})
    for up_document in uniprot_docs:
        evidence_level.append(up_document['data']['EVIDENCE']['LEVEL'])
    if len(evidence_level) > 0:
        median_evidence = np.median(evidence_level)
    else:
        median_evidence = np.nan

    data['nACCs'].append(len(uniref_dt['ACC']))
    data['nUniRef100'].append(len(uniref_dt['UNIREF']['UniRef100']))
    data['nUniRef90'].append(len(uniref_dt['UNIREF']['UniRef90']))
    data['nAF2'].append(len(uniref_dt['DARKNESS']['pLDDTs']))
    data['unirefID'].append(uniref_id)
#     data['uniref_rep_len'].append(rep_length)
    data['AF2_longest_best70'].append(longest_pLDDT_best70)
    data['AF2_longest_best70_len'].append(best_n_res)
    data['AF2_longest_best70_pLDDT'].append(best_pLDDT)
    data['median_Evidence'].append(median_evidence)

    if curr_count % write_step == 0:
        data = write_data_to_file(data, target_uniref, outfile)

        numb_seconds = time.time() - start
        time_to_end = round(((numb_seconds/curr_count)*n_entries)-numb_seconds)

        print('{} out of {}'.format(curr_count, n_entries), 'Time passed since start: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))), 'Expected to finish in: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(time_to_end)))-1, int(time.strftime('%d', time.gmtime(time_to_end)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(time_to_end))))


data = write_data_to_file(data, target_uniref, outfile)
            
numb_seconds = time.time() - start
print('\nFINISHED AFTER: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

