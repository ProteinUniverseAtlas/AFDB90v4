path_to_dbuilder = #path to dbuilder

import sys
import os

sys.path.append(path_to_dbuilder)
from src import extract_uniref    as uniref
from src import extract_interpro  as interpro
from src import extract_uniparc   as uniparc
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

interpro_db = interpro.interpro_db_diggested(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
interpro_db.index_db()

uniparc_db = uniparc.uniparc_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniparc_db.index_db()

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniprot_db.index_db()

# COUNT ENTRIES IN THE UNIREF DATABASE

print('\nCOUNTING DATABASE SIZE')

n_entries = uniref_db.col.count_documents({})
print(' ... FOUND {} ENTRIES'.format(n_entries))

# DEFINE OUTPUT FILE AND CHECK THE UNIREF IDS ALREADY THERE

print('\nDEFINING OUTPUT FILE')

outfile = 'AFDBv4_DUF_dark_diggestion_{}_2023-02-06.csv'.format(target_uniref)

# GO THROUGH EACH CHUNCK AND COLLECT THE TARGET DATA

start = time.time()

print('\nGOING THROUGH THE CHUNCKS')

write_step = 50000

data = {'unirefID': [], 'Has_duf': []}

curr_count = 0
for document in uniref_db.col.find():
    uniref_id = document['_id']
    uniref_dt = document['data']
    
    curr_count += 1
    
    if uniref_dt['DARKNESS']['FULL_noDUF'] <= 5:
        rep = uniref_dt['DARKNESS']['REP']
        has_duf = 0
        domains = []
        
        if rep is not None:
            if not rep.startswith('UP'):
                try:
                    domains = interpro_db.query(rep)[0]['data']
                except:
                    pass
                
                uniprot_dt = uniprot_db.query(rep)[0]['data']
                if 'CHAINS' in uniprot_dt:
                    domains += uniprot_dt['CHAINS']
                    
            else:
                try:
                    domains = uniparc_db.query(rep)[0]['data']['ANNO']
                except:
                    pass
                
            if len(domains) > 0:
                for domain in domains:
                    if 'DUF' in domain[0]:
                        has_duf = 1
            
        data['unirefID'].append(uniref_id)
        data['Has_duf'].append(has_duf)

    if curr_count % write_step == 0:
        data = write_data_to_file(data, target_uniref, outfile)

        numb_seconds = time.time() - start
        time_to_end = round(((numb_seconds/curr_count)*n_entries)-numb_seconds)

        print('{} out of {}'.format(curr_count, n_entries), 'Time passed since start: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))), 'Expected to finish in: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(time_to_end)))-1, int(time.strftime('%d', time.gmtime(time_to_end)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(time_to_end))))


data = write_data_to_file(data, target_uniref, outfile)
            
numb_seconds = time.time() - start
print('\nFINISHED AFTER: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))

