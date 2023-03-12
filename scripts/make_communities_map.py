import sys
import os
import json
import networkx as nx
import numpy as np
import pandas as pd
import time

infolder = sys.argv[1]

ingraph = '{}/full_graph.gml'.format(infolder)
node_data = '{}/node_class.json'.format(infolder)

def get_node_index(data, target_level = 'communityID', min_size = ):
    
    node_index = {}
    
    for i, node in enumerate(data['node']):
        node_index[node] = data[target_level][i]
    
    return node_index

def colapse_graph(graph, node_index):
    
    new_nodes = set(node_index.values())
    new_edges = set()
    
    n_expected = len(graph.edges)
    
    start = time.time()
    
    count = 0
    for edge in graph.edges:
        new_edge = [node_index[edge[0]], node_index[edge[1]]]
        if len(set(new_edge)) > 1:
            new_edge = tuple(sorted(new_edge))
            new_edges.add(new_edge)
        
        count+=1
        
        if count % 1000000 == 0:
            numb_seconds = time.time() - start
            time_to_end = round(((numb_seconds/count)*n_expected)-numb_seconds)
            print(count, n_expected, ' ... Time passed: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))), 'Expected to finish in: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(time_to_end)))-1, int(time.strftime('%d', time.gmtime(time_to_end)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(time_to_end))), flush = True)  
    
    print(' ... Defined {} edges connecting {} nodes'.format(len(new_edges), len(new_nodes)))
    
    return new_edges

def write_edges_list(edges, infolder):
    
    print('Writing edges file for cosmograph')
    
    start = time.time()
    
    outfile = '{}/communities_edge_list.csv'.format(infolder)
    
    with open(outfile, 'w') as outp:
        outp.write('innode,outnode\n')
        for edge in edges:
            if edge[0] != edge[1]:
                outp.write('{},{}\n'.format(edge[0], edge[1]))
  
    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
  
    
# load edges attributes
data = json.load(open(node_data, 'r'))

# load graph
graph = nx.read_gml(ingraph)

# get node index
node_index = get_node_index(data)

# colapse nodes
new_edges = colapse_graph(graph, node_index)

# write edges file for cosmograph
write_edges_list(new_edges, infolder)