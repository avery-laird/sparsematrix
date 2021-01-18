import pytest
from opt_ast import *

from scipy.io import mmread
from scipy.sparse import find, tril
import numpy
import subprocess
import networkx as nx
import matplotlib.pyplot as plt

from symbolicanalysis import reachset, level_order_set
from transformation import Prune

from utils import reachFromFiles, make_naive_solve

from gen_triang import generate_c_parallel


def test_reachset_small():
    r, _, _ = reachFromFiles('rset_example.mtx', 'rset_example_b.mtx')
    assert list(nx.topological_sort(r)) == [6, 1, 7, 8, 9, 10] or list(nx.topological_sort(r) == [1, 6, 7, 8, 9, 10])


# def test_reachset_torso():
#     r = reachFromFiles('torso1/torso1.mtx', 'b_for_torso1.mtx')
#     assert True


def test_int_decl():
    int_decl_p = Declaration(IntType(), Iden("p"))
    assert int_decl_p.gen() == "int p;"

    int_decl_j = Declaration(IntType(), Iden("j"))
    assert int_decl_j.gen() == "int j;"


def test_assign():
    j_sym = Iden("j")
    rhs = Int(0)
    assign_j_rhs = Assign(j_sym, rhs)
    assert assign_j_rhs.gen() == "j = 0"


def test_cmp_less():
    j_sym = Iden("j")
    n_sym = Iden("n")
    less = Less(j_sym, n_sym)
    assert less.gen() == "j < n"


def test_address():
    iden1 = Iden("Lx")
    iden2 = Iden("Lp")
    addr2 = Address(iden2, Iden("j"))
    addr1 = Address(iden1, addr2)
    assert addr1.gen() == "Lx[Lp[j]]"


def test_postfix_inc():
    j_iden = Iden("j")
    assert PostfixInc(j_iden).gen() == "j++"


def test_forloop():
    j_iden = Iden("j")
    zero = Int(0)
    init = Assign(j_iden, zero)
    n_iden = Iden("n")
    cond = Less(j_iden, n_iden)
    do = PostfixInc(j_iden)
    forloop = ForLoop([init, cond, do], [])
    assert forloop.gen() == """for (j = 0; j < n; j++) {

}"""


def test_naive_solve():
    top, _ = make_naive_solve()
    assert top.gen() == """for (j = 0; j < n; j++) {
x[j] /= Lx[Lp[j]];
for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
x[Li[p]] -= Lx[p] * x[j];
}
}"""


def test_reachset_transform():
    rset, _, _ = reachFromFiles('rset_example.mtx', 'rset_example_b.mtx')
    outer, inner = make_naive_solve()
    reachSetInit, loop = Prune(outer, rset).run()

    assert loop.gen() == """for (px = 0; px < 6; px++) {
j = reachSet[px];
x[j] /= Lx[Lp[j]];
for (p = Lp[j] + 1; p < Lp[j + 1]; p++) {
x[Li[p]] -= Lx[p] * x[j];
}
}"""


def test_reachset_torso1():
    rset, _, _ = reachFromFiles('torso1/torso1.mtx', 'b_for_torso1.mtx')
    outer, inner = make_naive_solve()
    reachSetInit, loop = Prune(outer, rset).run()
    assert True


def test_level_order_set():
    rset, A, b = reachFromFiles('rset_example.mtx', 'rset_example_b.mtx')
    level_set, etree = level_order_set(A)
    assert level_set == [
        [0, 1, 3],
        [6, 2, 4],
        [5],
        [7],
        [8],
        [9]
    ]
    test = nx.DiGraph()
    test.add_edges_from([
        (0, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (5, 7),
        (4, 5),
        (3, 4),
        (2, 5),
        (1, 2)
    ])
    assert nx.is_isomorphic(etree, test)


def test_generate_c_parallel():
    output = generate_c_parallel('rset_example.mtx', 'rset_example_b.mtx', '2')
    print(output)
