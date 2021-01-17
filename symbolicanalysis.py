from __future__ import annotations

from scipy.sparse.coo import coo_matrix
from scipy.sparse import tril

import networkx as nx
import matplotlib.pyplot as plt


# class Node:
#     def __init__(self, label):
#         self.label = label
#         self.outgoing: list[Node] = []
#
#     def to(self, other: Node):
#         self.outgoing.append(Node)
#
# class Graph:
#     def __init__(self):
#         self.nodes: list[Node] = []


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
    A = matrix
    (m, n) = A.shape
    G = nx.DiGraph()
    G.add_nodes_from(range(1, n + 1))
    (rows, cols) = tril(A, -1).nonzero()
    for i, j in zip(rows, cols):
        G.add_edge(j + 1, i + 1)
    # do a DFS from beta
    reach = nx.DiGraph()
    beta = b
    (betaRows, _) = beta.nonzero()
    for r in betaRows:
        tree = nx.dfs_tree(G, r + 1)
        reach = nx.compose(reach, tree)

    # nx.draw(reach, with_labels=True)
    # plt.show()
    return reach
