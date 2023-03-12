import sys
import os
import json
import networkx as nx
import numpy as np
import time
import gzip
import io

from networkx.algorithms.community import asyn_lpa_communities

# LOAD INPUTS

infasta = sys.argv[1]
inmmsqs = sys.argv[2]
outfolder = sys.argv[3]

if not os.path.isdir(outfolder):
    os.mkdir(outfolder)

# helping routines

def get_seqs_index(infasta, indx = {}):

    print('\nReading sequences from input fasta and generating node index')
    
    start = time.time()
    
    count = len(indx)
    
    if infasta.endswith('.gz'):
        with gzip.open(infasta, 'rb') as inf:
            with io.TextIOWrapper(inf, encoding='utf-8') as decoder:
                for line in decoder:
                    if line.startswith('>'):
                        if '|' in line:
                            line = line.split('|')[1].strip('>')
                        else:
                            line = line.split()[0].strip('>')
                        indx[line] = count
                        count+=1
    
    else:
        with open(infasta, 'r') as inf:
            for line in inf:
                if line.startswith('>'):
                    if '|' in line:
                        line = line.split('|')[1].strip('>')
                    else:
                        line = line.split()[0].strip('>')
                    indx[line] = count
                    count+=1

    print(' ... No. of expected nodes:', len(indx))
    
    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
  
    return indx


def get_neighbors_simple(inmmsqs, indexes, mineval = 1e-4, mincov = 0, simplex = True, nmax=None):
    
    print('\nCollecting edges from input mmseqs file')
    
    start = time.time()
    
    edges = set()
    nodes = set()
    
    if not simplex:
        edges = {}
    
    mincov = mincov/100
    previous_len = None
    with open(inmmsqs, 'r') as inmm:
        for line in inmm:
            line = line.split('\t')
            i, j, cov, evalue = line[0], line[1], line[2], line[10]

            if i in indexes and j in indexes:
                evalue = float(evalue.strip())
                cov = float(cov)

                if i != j and evalue <= mineval and cov >= mincov: 
                    if simplex:
#                         edges.add(tuple(sorted([indexes[i], indexes[j]])))
                        edges.add(tuple(sorted([i, j])))
                    else:
#                         edge = sorted([indexes[i], indexes[j]])
                        edge = sorted([i, j])
                        if edge[0] in edges:
                            if edge[1] in edges[edge[0]]:
                                if edges[edge[0]][edge[1]][0] > evalue:
                                    edges[edge[0]][edge[1]] = (evalue, cov)
                            else:
                                edges[edge[0]][edge[1]] = (evalue, cov)
                        else:
                            edges[edge[0]] = {edge[1]: (evalue, cov)}
                            
                    nodes.add(i)
                    nodes.add(j)
            
                if len(nodes)> 0 and len(nodes) % 100000 == 0 and len(nodes) != previous_len:
                    print(len(nodes))
                    previous_len = len(nodes)
                
                if nmax is not None and len(nodes) == nmax:
                    break
    
    print(' ... Total number of hubs:',  len(edges))
    print(' ... Total number of nodes:', len(nodes))
    
    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
    
    return edges, nodes

def generate_pairs(neigbhrs, topN = None, min_weight = 0):
    
    print('\nRemoving redundant edges')
    
    start = time.time()
    
    edges = dict()
    weights = list()
    
    for i in neigbhrs:
        if topN is not None:
            curr_neighbrs = {k: v for k, v in sorted(neigbhrs[i].items(), key=lambda item: item[1])[:topN]}
        else:   
            curr_neighbrs = neigbhrs[i]
         
        for j in curr_neighbrs:
            i_index, j_index = sorted([i, j])
            evalue = curr_neighbrs[j][0]
            cov = curr_neighbrs[j][1]
            
            edge = (i_index, j_index)
            edges[edge] = {'evalue': evalue, 'cov': cov*100}
#             edges.add(edge)
#             weights.append(-np.log10(evalue)*cov)
    
    print(' ... Total number of edges:', len(edges))
    
#     # normalise weights
#     max_weigth = max(weights) 
#     min_weigth = int(min(weights))
#     normalised_weigths = [(weight-min_weight)/(max_weigth-min_weight) for weight in weights]
    
#     for i, edge in enumerate(edges):
#         edges[edge]['weight'] = normalised_weigths[i]
        
    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
  
    return edges

def build_graph(edges, indexes, outgraph = None, outfolder = outfolder, map_properties = False, properties=['Darkness']):
    
    print('\nBuilding the graph')
    
    start = time.time()
    
    G=nx.Graph()
    G.add_nodes_from(list(indexes.keys()))
    G.add_edges_from(list(edges.keys()))
    
    nx.set_edge_attributes(G, edges)
    
    if map_properties:
        properties = get_nodes_properties(indexes.keys(), properties)
        nx.set_node_attributes(G, properties)
    
    if outgraph is None:
        nx.write_gml(G, "{}/full_graph.gml".format(outfolder))
    else:
        nx.write_gml(G, outgraph)
    
    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
    
    return G
    
    
def collect_connected_components(G, nodes, min_size = 0, outfolder=outfolder, outgraph=None):
    
    print('\nCollecting individual subgraphs/connected components')
    
    sec_outfolder = '{}/subgraphs'.format(outfolder)
    if not os.path.isdir(sec_outfolder):
        os.mkdir(sec_outfolder)
    
    start = time.time()
            
    # Get connected components
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    
    print(' ... Found {} subgraphs'.format(nx.number_connected_components(G)))
    print(' ... ... The largest has {} nodes'.format(len(components[0])))
    
    count = 0
    edge_count = 0
    node_count = 0
    
    node_cluster_class = {'node': [], 'subgraphID': [], 'communityID': []}

    for component_index, c in enumerate(components):
    
        curr_size = len(c)

        if curr_size >= min_size:
            
            component = G.subgraph(c).copy()
            
            curr_outf = '{}/subgraph_{:06d}.gml'.format(sec_outfolder, component_index)
            nx.write_gml(component, curr_outf)

#             Get communities by label propagation
            communities = list(asyn_lpa_communities(component, weight=None))

            for community_index, community in enumerate(communities):
                node_cluster_class['node'] += community
                node_cluster_class['communityID'] += ['{}[{}]'.format(component_index, community_index) for node in community]
                node_cluster_class['subgraphID'] += [component_index for node in community]
            
            print(' ... ... Subgraph:', component_index, 'No. nodes:', curr_size, 'No. communities:', len(communities))
            print(' ... ... Subgraph:', component_index, 'No. nodes:', curr_size)

            
    json.dump(node_cluster_class, open('{}/node_class.json'.format(outfolder), 'w'))
    
#     if outgraph is None:
#         nx.write_gml(G, "{}/full_graph.gml".format(outfolder))
#     else:
#         nx.write_gml(G, outgraph)
    
    print(' ... Wrote {} subgraphs, totalling {} nodes and {} edges'.format(count, node_count, edge_count))
    
    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
    

# MAIN CODE

outgraph = "{}/full_graph.gml".format(outfolder)

indexes = get_seqs_index(infasta)

if not os.path.isfile(outgraph):
    hubs, nodes = get_neighbors_simple(inmmsqs, indexes, mineval = 1e-4, mincov = 50, simplex = False, nmax=None)
    edges = generate_pairs(hubs, topN = 4)

    graph = build_graph(edges, indexes, outgraph=outgraph)

else:
    print('Graph already produced. Will just load it')
    
    start = time.time()
    graph = nx.read_gml(outgraph)

    numb_seconds = time.time() - start
    print(' ... Took me: {} months {} days {}'.format(int(time.strftime('%m', time.gmtime(numb_seconds)))-1, int(time.strftime('%d', time.gmtime(numb_seconds)))-1, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))


collect_connected_components(graph, nodes=list(indexes.keys()), min_size = 2)
