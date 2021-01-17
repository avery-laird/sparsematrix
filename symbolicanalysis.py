from __future__ import annotations

from scipy.sparse.coo import coo_matrix
from scipy.sparse import tril

import networkx as nx
import matplotlib.pyplot as plt


def initial_directed_graph(matrix: coo_matrix) -> nx.DiGraph:
    """
    Generate a DAG from matrix A
    :param matrix:
    :param b:
    :return:
    """
    (m, n) = matrix.shape
    G = nx.DiGraph()
    G.add_nodes_from(range(1, n + 1))
    (rows, cols) = tril(matrix, -1).nonzero()
    for i, j in zip(rows, cols):
        G.add_edge(j + 1, i + 1)
    return G


def reachset(matrix: coo_matrix, b: coo_matrix) -> nx.DiGraph:
    """
    computes the reach set, the set of all nodes
    reachable from any node in beta = {i | b_i != 0}
    First obtain DG_L:
        create nodes 1 .. n, n = matrix rank, but it seems # col is usually used
        create edge (j,i) if L[i][j] != 0
    Then search b, and put each non-zero row index in beta
    :param matrix:
    :param b:
    :return:
    """
    G = initial_directed_graph(matrix)
    # do a DFS from beta
    reach = nx.DiGraph()
    beta = b
    (betaRows, _) = beta.nonzero()
    for r in betaRows:
        tree = nx.dfs_tree(G, r + 1)
        reach = nx.compose(reach, tree)
    return reach


def level_order_set(A: coo_matrix):
    """
    Given a matrix A, return its elimination tree
    :param reachset:
    :param A:
    :return:
    """
    csc = tril(A).tocsc()
    etree = nx.DiGraph()
    # add an edge for each first non-zero value
    # in the row
    nonzeros = len(csc.indices)
    for col, col_start in enumerate(csc.indptr):
        if col_start + 1 >= nonzeros: break;
        # first off-diagonal non-zero entry
        etree.add_edge(col, csc.indices[col_start + 1])
    # for each source in the etree, do a BFS
    level_set = []
    schedule = etree.copy()
    while schedule.order() > 0:
        level = []
        for node in list(schedule.nodes):
            if schedule.in_degree(node) > 0: continue
            # otherwise, add the node to the level set
            # and remove it from the etree
            level.append(node)
        schedule.remove_nodes_from(level)
        level_set.append(level)
    return level_set, etree