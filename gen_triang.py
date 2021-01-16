from transformation import Prune
from opt_ast import *
from utils import reachFromFiles, make_naive_solve
import networkx as nx

optsolve = """int optsolve(int n, int *Lp, int *Li, double *Lx, double *&x) {{
    int px, p, j;
    if (!Lp || !Li || !x) return (0); /* check inputs */
    {}
    return (1);
}}
"""

def generate_c(apath, bpath):
    rset: nx.DiGraph = reachFromFiles(apath, bpath)
    outer, inner = make_naive_solve()
    reachSetInit, loop = Prune(outer, rset).run()
    return optsolve.format("\n".join([reachSetInit.gen(), loop.gen()]))

if __name__ == "__main__":
    print(generate_c('rset_example.mtx', 'rset_example_b.mtx'))
    # print(generate_c('torso1/torso1.mtx', 'b_for_torso1.mtx'))