from transformation import Prune
from opt_ast import *
from utils import reachFromFiles, make_naive_solve
import networkx as nx
import sys

optsolve = """int optsolve{}(int n, int *Lp, int *Li, double *Lx, double *&x) {{
    int px, p, j;
    if (!Lp || !Li || !x) return (0); /* check inputs */
    {}
    return (1);
}}
"""


def generate_c(apath, bpath, funcid):
    rset: nx.DiGraph = reachFromFiles(apath, bpath)
    outer, inner = make_naive_solve()
    reachSetInit, loop = Prune(outer, rset).run()
    return optsolve.format(funcid, "\n".join([reachSetInit.gen(), loop.gen()]))


def test_helper(paths):
    for i, (A, b) in enumerate(paths):
        print("Generating test {} ({}) into optimized_triang_test{}.cpp...".format(i, A, i))
        with open("optimized_triang_test{}.cpp".format(i), 'w') as f:
            f.write(generate_c(A, b, str(i)))


if __name__ == "__main__":
    if (len(sys.argv) == 2):
        # then use the arguments to generate a new file
        # the expected format is python gen_triang.py <A path> <b path>
        print("Generating test from {} and printing to stdout (redirect to whatever file is convenient)".format(sys.argv[1]))
        print(generate_c(sys.argv[1], sys.argv[2], "_"))
        exit(0)

    test_helper([
        # ('torso1/torso1.mtx', 'b_for_torso1.mtx'),
        # ('TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx', 'b_for_TSOPF_RS_b678_c2_b.mtx'),
        ('rset_example.mtx', 'rset_example_b.mtx')
    ])

    # print("Generating test 1 (torso1.mtx) into optimized_triang_test1.cpp...")
    # with open("optimized_triang_test1.cpp") as f:
    #     f.write(generate_c('torso1/torso1.mtx', 'b_for_torso1.mtx'))
    #
    # print("Generating test 2 (TSOPF_RS_b678_c2.mtx) into optimized_triang_test2.cpp...")
    # with open("optimized_triang_test2.cpp") as f:
    #     f.write(generate_c('TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx', 'b_for_TSOPF_RS_b678_c2_b.mtx'))
    #
    # print("Generating test 3 (TSOPF_RS_b678_c2.mtx) into optimized_triang_test2.cpp...")
    # with open("optimized_triang_test2.cpp") as f:
    #     f.write(generate_c('TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx', 'b_for_TSOPF_RS_b678_c2_b.mtx'))
    # print(generate_c('rset_example.mtx', 'rset_example_b.mtx'))
