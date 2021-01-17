import pytest
from opt_ast import *

from scipy.io import mmread
from scipy.sparse import find, tril
import numpy
import subprocess
import networkx as nx

from symbolicanalysis import reachset
from transformation import Prune

from utils import reachFromFiles, make_naive_solve


# class Solving(unittest.TestCase):
#     def setUp(self):
#         # build starter.cpp
#         subprocess.run(["cmake", "-B", "cmake-build-debug"], check=True)
#         subprocess.run(["make"],
#                        cwd="/home/avery/Projects/sparse_matrix_opt/cmake-build-debug", check=True)
#
#     def test_correct_result(self):
#         # run starter.cpp
#         # subprocess.run([
#         #     "/home/avery/Projects/sparse_matrix_opt/cmake-build-debug/starter"
#         # ])
#         A = tril(mmread('torso1/torso1.mtx'))
#         b = mmread('b_for_torso1.mtx')
#         x = mmread('testx.mtx')
#         xcheck = A.multiply(x)
#         # https://stackoverflow.com/questions/47770906/how-to-test-if-two-sparse-arrays-are-almost-equal
#         (i1, j1, k1) = find(xcheck)
#         (i2, j2, k2) = find(b)
#         return \
#             numpy.array_equal(i1, i2) and \
#             numpy.array_equal(j1, j2) and \
#             numpy.array_equal(k1, k2)


def test_reachset_small():
    r = reachFromFiles('rset_example.mtx', 'rset_example_b.mtx')
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
    rset: nx.DiGraph = reachFromFiles('rset_example.mtx', 'rset_example_b.mtx')
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
    rset: nx.DiGraph = reachFromFiles('torso1/torso1.mtx', 'b_for_torso1.mtx')
    outer, inner = make_naive_solve()
    reachSetInit, loop = Prune(outer, rset).run()

    assert True
