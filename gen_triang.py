from transformation import Prune, Parallelize
from opt_ast import *
from utils import reachFromFiles, make_naive_solve, levelOrderFromFiles
import networkx as nx
import sys

optsolve = """int optsolve{}(int n, int *Lp, int *Li, double *Lx, double *&x) {{
    int px, p, j;
    if (!Lp || !Li || !x) return (0); /* check inputs */
    {}
    return (1);
}}
"""

# contains 5 points of insertion:
# 1. function name suffix
# 2. level set initialization
# 3. level set length
# 4. lengths of level set elements
# 5. parallel for loop
optsolve_par = """#include <omp.h>
#include <vector>

int optsolve{}(int n, int *Lp, int *Li, double *Lx, double *&x) {{
    int px, p, j;
    /* level set, level_sets */
    {}
    /* length, level_set_len */
    int level_set_len;
    {}
    /* length array, level_set_element_length */
    {}
    if (!Lp || !Li || !x) return (0); /* check inputs */
    for (int i=0; i<level_set_len; ++i) {{
#pragma omp parallel for
    {}
}} /* end outer for loop */
    return (1);
}}
"""


def generate_c(apath, bpath, funcid):
    rset, A, b = reachFromFiles(apath, bpath)
    outer, inner = make_naive_solve()
    reachSetInit, loop = Prune(outer, rset).run()
    return optsolve.format(funcid, "\n".join([reachSetInit.gen(), loop.gen()]))


def generate_c_parallel(apath, bpath, funcid):
    level_order, _ = levelOrderFromFiles(apath, bpath)
    outer, inner = make_naive_solve()
    transformed = Parallelize(outer, level_order).run()
    return optsolve_par.format(
        funcid + "_parallel",
        transformed.level_set_init.gen(),
        transformed.level_set_length.gen(),
        transformed.level_set_element_length_init.gen(),
        transformed.for_loop.gen(),
    )


def test_helper(paths, parallel=False):
    generator = generate_c_parallel if parallel else generate_c
    suffix = "_parallel" if parallel else ""
    for i, (A, b) in enumerate(paths):
        print("Generating test {} ({}) into optimized_triang_test{}.cpp...".format(i, A, str(i) + suffix))
        with open("optimized_triang_test{}.cpp".format(str(i) + suffix), 'w') as f:
            f.write(generator(A, b, str(i)))


if __name__ == "__main__":
    if (len(sys.argv) == 3):
        # then use the arguments to generate a new file
        # the expected format is python gen_triang.py <A path> <b path>
        print("Generating test from {} and printing to stdout (redirect to whatever file is convenient)".format(
            sys.argv[1]))
        print(generate_c(sys.argv[1], sys.argv[2], "_"))
        exit(0)

    test_helper([
        ('torso1/torso1.mtx', 'b_for_torso1.mtx'),
        ('TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx', 'b_for_TSOPF_RS_b678_c2_b.mtx'),
        ('rset_example.mtx', 'rset_example_b.mtx')
    ])

    test_helper([
        ('torso1/torso1.mtx', 'b_for_torso1.mtx'),
        ('TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx', 'b_for_TSOPF_RS_b678_c2_b.mtx'),
        ('rset_example.mtx', 'rset_example_b.mtx')
    ], parallel=True)
