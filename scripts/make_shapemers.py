import argparse
import gzip
import io
import itertools

import prody as pd
import tarfile
import torch
from scipy import ndimage
from time import time

from geometricus import MultipleMomentInvariants, ShapemerLearn

from pathlib import Path
import json
import numba as nb
import numpy as np
from tqdm import tqdm

import proteinnet_parser


def parse_pae_file(pae_json_data):
    if type(pae_json_data) == str or type(pae_json_data) == Path:
        with open(pae_json_data, "rt") as f:
            data = json.load(f)[0]
    else:
        data = json.load(pae_json_data)[0]

    if 'residue1' in data and 'distance' in data:
        # Legacy PAE format, keep for backwards compatibility.
        r1, d = data['residue1'], data['distance']
        size = max(r1)
        matrix = np.empty((size, size), dtype=np.float64)
        matrix.ravel()[:] = d
    elif 'predicted_aligned_error' in data:
        # New PAE format.
        matrix = np.array(data['predicted_aligned_error'], dtype=np.float64)
    else:
        raise ValueError('Invalid PAE JSON format.')

    return matrix


@nb.njit
def get_plddt_matrix(plddt):
    size = len(plddt)
    matrix = np.empty((size, size), dtype=np.float64)
    for i in range(size):
        matrix[i, i] = plddt[i]
        for j in range(i + 1, size):
            matrix[i, j] = matrix[j, i] = (plddt[i] + plddt[j])
    return 100 - matrix / 2


def get_domains_networkx(pae_matrix, plddt_matrix, cutoff=20, graph_resolution=0.5):
    """
    Adapted from https://github.com/tristanic/pae_to_domains
    
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero value
        to avoid divide-by-zero warnings
        * plddt_matrix: a (n_residues x n_residues) numpy array containing average pairwise PLDDT values
        * cutoff (optional, default=20): graph edges will only be created for residue pairs with pae+avg.PLDDT <
        cutoff
        * graph_resolution (optional, default=0.5): regulates how aggressively the clustering algorithm is.
        Smaller values lead to larger clusters. Value should be larger than zero, and values larger than 5 are
        unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    """
    try:
        import networkx as nx
        from networkx.algorithms import community
    except ImportError:
        print(
            'ERROR: This method requires NetworkX (>=2.6.2) to be installed. Please install it using "pip install '
            'networkx" in a Python >=3.7 environment and try again.')
        import sys
        sys.exit()
    matrix = pae_matrix + plddt_matrix
    weights = 1 / matrix
    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    edges = np.argwhere(matrix < cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    wedges = [(i, j, w) for (i, j), w in zip(edges, sel_weights)]
    g.add_weighted_edges_from(wedges)
    clusters = community.greedy_modularity_communities(g, weight='weight', resolution=graph_resolution)
    return clusters


def get_domains_igraph(pae_matrix, plddt_matrix, cutoff=20, graph_resolution=0.5):
    """
    Adapted from https://github.com/tristanic/pae_to_domains
    
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * plddt_matrix: a (n_residues x n_residues) numpy array containing average pairwise PLDDT values
        * cutoff (optional, default=20): graph edges will only be created for residue pairs with pae+avg.PLDDT<cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    """
    try:
        import igraph
    except ImportError:
        print(
            'ERROR: This method requires python-igraph to be installed. Please install it using "pip install '
            'python-igraph" in a Python >=3.6 environment and try again.')
        import sys
        sys.exit()
    matrix = pae_matrix + plddt_matrix
    weights = 1 / matrix
    g = igraph.Graph()
    size = weights.shape[0]
    g.add_vertices(range(size))
    edges = np.argwhere(pae_matrix < cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    g.add_edges(edges)
    g.es['weight'] = sel_weights

    vc = g.community_leiden(weights='weight', resolution_parameter=graph_resolution / 100, n_iterations=-1)
    membership = np.array(vc.membership)
    from collections import defaultdict
    clusters = defaultdict(list)
    for i, c in enumerate(membership):
        clusters[c].append(i)
    return clusters.values()


def clusters_to_domains(protein, clusters, min_length=20, avg_plddt_cutoff=70):
    chain = "A"
    for cl in clusters:
        start, stop = min(cl), max(cl)
        if (stop - start) >= min_length:
            domain = protein.select(f"resnum {start}:{stop}")
            if domain.getBetas().mean() >= avg_plddt_cutoff:
                for res in domain:
                    res.setChid(chain)
                yield domain
                chain = chr(ord(chain) + 1)


def split_alphafold_protein(prody_protein, pae_file=None, plddt_threshold=70, sigma=5):
    """
    Splits an AlphaFold protein into fragments based on a Gaussian-smoothed version of the PLDDT score.
    Parameters
    ----------
    prody_protein
        ProDy protein object of calpha atoms
    pae_file
        pae_file_data
    plddt_threshold
        Fragments will be split according to residues with a (smoothed) PLDDT score below this threshold.
    sigma
        Sigma for the smoothing of the PLDDT score.

    Returns
    -------
    (start, end) indices for each split
    """
    if pae_file is not None:
        pae_matrix = parse_pae_file(pae_file)
        beta_list = prody_protein.getBetas()
        plddt_matrix = get_plddt_matrix(beta_list)
        clusters = get_domains_igraph(pae_matrix, plddt_matrix)
        beta_list = ndimage.gaussian_filter1d(beta_list, sigma=sigma)
        all_slices = []
        for cl in clusters:
            chain_start, chain_stop = min(cl), max(cl)
            length = chain_stop - chain_start
            if length < 20:
                continue
            indices = np.ones(length, dtype=int)
            indices[np.where(beta_list[chain_start:chain_stop] < plddt_threshold)] = 0
            slices = ndimage.find_objects(ndimage.label(indices)[0])
            slices = [(s[0].start, s[0].stop) for s in slices]
            all_slices += [(chain_start + start, chain_start + stop) for start, stop in slices]
    else:
        beta_list = ndimage.gaussian_filter1d(prody_protein.getBetas(), sigma=sigma)
        indices = np.ones(beta_list.shape[0], dtype=int)
        indices[np.where(beta_list < plddt_threshold)] = 0
        slices = ndimage.find_objects(ndimage.label(indices)[0])
        all_slices = [(s[0].start, s[0].stop) for s in slices]
    return all_slices


def split_pdb_protein(prody_protein):
    """
    Splits a protein into fragments based on chain.
    Parameters
    ----------
    prody_protein
        ProDy protein object.

    Returns
    -------
    (start, end, chid) indices for each split
    """
    slices = []
    chains = set(a.getChid() for a in prody_protein)
    if len(chains):
        for chain in chains:
            if not len(chain.strip()):
                chain = prody_protein
            else:
                chain = prody_protein.select(f"chain {chain}")
            slices.append((chain[0].getResindex(), chain[-1].getResindex() + 1, chain[0].getChid()))
    else:
        slices.append((prody_protein[0].getResindex(), prody_protein[-1].getResindex(), ''))
    return sorted(slices)


def get_shapemers(calpha_protein,
                  model,
                  is_af=False,
                  pae_file_data=None,
                  length_threshold=20,
                  plddt_threshold=70,
                  sigma=5):
    """
    Retrieves the moments of the protein.
    Parameters
    ----------
    calpha_protein
        prody object
    model
        ShapemerLearn model
    is_af
        Whether the protein is an AlphaFold protein
    pae_file_data
        pae file as extrected gzip or as filename
    length_threshold
        Proteins with fewer (filtered) residues than this threshold will be ignored.
    plddt_threshold
        Residues with a (smoothed) PLDDT score below this threshold will be ignored.
    sigma
        Sigma for the smoothing of the PLDDT score.

    Returns
    -------
    """
    if is_af:
        residue_slices = split_alphafold_protein(calpha_protein, pae_file_data, plddt_threshold, sigma)
    else:
        residue_slices = split_pdb_protein(calpha_protein)
    coords = calpha_protein.getCoords()
    return get_shapemers_from_coords(coords, model, length_threshold=length_threshold, residue_slices=residue_slices)


def get_shapemers_from_coords(coords, model, length_threshold=20, residue_slices=None):
    shapemers = []
    indices = [len(coords)]
    if residue_slices is None:
        residue_slices = [(0, len(coords))]
    try:
        for x in residue_slices:
            start_index, end_index, *_ = x
            if end_index - start_index > length_threshold:
                indices += list(range(start_index, end_index))
                shapemers += MultipleMomentInvariants.from_coordinates("name",
                                                                       coords[
                                                                       start_index:end_index]).get_shapemers_model(
                    model)
        if len(shapemers):
            assert len(shapemers) == len(indices) - 1
            return indices, shapemers
    except Exception as e:
        print(f"Error {e}")
    return [], []


def make_corpus_proteome(taxid, db_folder, output_folder):
    model = ShapemerLearn.load()
    shapemer_keys = list(map(tuple, itertools.product([0, 1], repeat=model.output_dimension)))
    key_to_index = dict(zip(shapemer_keys, range(len(shapemer_keys))))
    start = time()
    index = 0
    f_s = open(output_folder / f"{taxid}_shapemers.txt", "w")
    f_i = open(output_folder / f"{taxid}_indices.txt", "w")
    for f in db_folder.glob(f"proteome-tax_id-{taxid}-*_v4.tar"):
        with tarfile.open(f) as tar:
            for fh in tar.getmembers():
                if '.cif' in fh.name:
                    if index % 1000 == 0:
                        print(f"{index} proteins processed in {time() - start} seconds")
                    uniprot_ac = '-'.join(fh.name.split('-')[1:3])
                    with io.TextIOWrapper(gzip.open(tar.extractfile(fh), 'r'), encoding='utf-8') as mmcif:
                        with gzip.open(tar.extractfile(
                                tar.getmember(f"AF-{uniprot_ac}-predicted_aligned_error_v4.json.gz"))) as pae:
                            protein = pd.parseMMCIFStream(mmcif)
                            protein = protein.select("protein and calpha")
                            indices, shapemers = get_shapemers(protein, model, is_af=True, pae_file_data=pae)
                            if len(shapemers):
                                f_i.write(f"{uniprot_ac}\t{indices[0]}\t{' '.join(str(s) for s in indices[1:])}\n")
                                f_s.write(f"{uniprot_ac}\t{' '.join(str(key_to_index[s]) for s in shapemers)}\n")
                    index += 1
    f_s.close()
    f_i.close()


def make_corpus_from_file(filename, db_folder, output_folder):
    model = ShapemerLearn.load()
    shapemer_keys = list(map(tuple, itertools.product([0, 1], repeat=model.output_dimension)))
    shapemer_key_to_index = dict(zip(shapemer_keys, range(len(shapemer_keys))))
    f_s = open(output_folder / f"{filename.stem}_shapemers.txt", "w")
    f_i = open(output_folder / f"{filename.stem}_indices.txt", "w")
    with open(filename) as f:
        num_lines = sum(1 for _ in f)
    with tarfile.open(db_folder / f"{filename.stem}.tar") as tar:
        with open(filename) as f:
            for line in tqdm(f, total=num_lines):
                fh = line.strip()
                uniprot_ac = '-'.join(fh.split('-')[1:3])
                with io.TextIOWrapper(gzip.open(tar.extractfile(tar.getmember(fh)), 'r'), encoding='utf-8') as mmcif:
                    with gzip.open(tar.extractfile(
                            tar.getmember(f"AF-{uniprot_ac}-predicted_aligned_error_v4.json.gz"))) as pae:
                        protein = pd.parseMMCIFStream(mmcif)
                        protein = protein.select("protein and calpha")
                        indices, shapemers = get_shapemers(protein, model, is_af=True, pae_file_data=pae)
                        if len(shapemers):
                            f_i.write(f"{uniprot_ac}\t{indices[0]}\t{' '.join(str(s) for s in indices[1:])}\n")
                            f_s.write(
                                f"{uniprot_ac}\t{' '.join(str(shapemer_key_to_index[s]) for s in shapemers)}\n")
    f_s.close()
    f_i.close()


def make_corpus_pdb_folder(db_folder_divided, output_folder):
    model = ShapemerLearn.load()
    shapemer_keys = list(map(tuple, itertools.product([0, 1], repeat=model.output_dimension)))
    key_to_index = dict(zip(shapemer_keys, range(len(shapemer_keys))))
    f_s = open(output_folder / f"{db_folder_divided.stem}_shapemers.txt", "w")
    f_i = open(output_folder / f"{db_folder_divided.stem}_indices.txt", "w")
    for filename in tqdm(db_folder_divided.glob(f"*.ent.gz")):
        try:
            uid = filename.stem.split(".")[0].split("pdb")[1]
            with gzip.open(filename, 'r') as pdb:
                with io.TextIOWrapper(pdb, encoding='utf-8') as decoder:
                    protein = pd.parsePDBStream(decoder)
            protein = protein.select("protein and calpha")
            if protein is None:
                continue
            coords = protein.getCoords()
            residue_slices = split_pdb_protein(protein)
            for start_index, end_index, chain in residue_slices:
                indices, shapemers = get_shapemers_from_coords(coords, model, residue_slices=[(start_index, end_index)])
                if len(shapemers):
                    f_i.write(f"{uid}_{chain}\t{indices[0]}\t{' '.join(str(s) for s in indices[1:])}\n")
                    f_s.write(f"{uid}_{chain}\t{' '.join(str(key_to_index[s]) for s in shapemers)}\n")
        except Exception as e:
            print(e)
    f_s.close()
    f_i.close()


def make_corpus_proteinnet(db_folder, output_folder):
    model = ShapemerLearn.load()
    shapemer_keys = list(map(tuple, itertools.product([0, 1], repeat=model.output_dimension)))
    key_to_index = dict(zip(shapemer_keys, range(len(shapemer_keys))))
    f_s = open(output_folder / f"proteinnet_shapemers.txt", "w")
    f_i = open(output_folder / f"proteinnet_indices.txt", "w")
    for filename in [db_folder / x for x in ["training_100", "validation", "testing"]]:
        with open(filename) as f:
            total = sum(1 for line in f if line == "[ID]\n")
        for entry in tqdm(proteinnet_parser.yield_records_from_file(filename, 20), total=total):
            entry = proteinnet_parser.clean_entry(entry, 'ca')
            uid = entry["ID"]
            indices, shapemers = get_shapemers_from_coords(entry["tertiary"], model)
            if len(shapemers):
                f_i.write(f"{uid}\t{indices[0]}\t{' '.join(str(s) for s in indices[1:])}\n")
                f_s.write(f"{uid}\t{' '.join(str(key_to_index[s]) for s in shapemers)}\n")
    f_s.close()
    f_i.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=Path)
    parser.add_argument("db_folder", type=Path)
    parser.add_argument("output_folder", type=Path)
    args = parser.parse_args()
    if args.filename:
        make_corpus_from_file(args.filename, args.db_folder, Path(str(args.output_folder).strip()))


def main_pdb():
    parser = argparse.ArgumentParser()
    parser.add_argument("db_folder_divided", type=Path)
    parser.add_argument("output_folder", type=Path)
    args = parser.parse_args()
    make_corpus_pdb_folder(args.db_folder_divided, Path(str(args.output_folder).strip()))


if __name__ == '__main__':
    main_pdb()
