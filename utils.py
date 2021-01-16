from scipy.io import mmread
from symbolicanalysis import reachset
from opt_ast import *


def make_naive_solve():
    x = Iden("x")
    Li = Iden("Li")
    p = Iden("p")
    Lx = Iden("Lx")
    j = Iden("j")
    Lp = Iden("Lp")
    n = Iden("n")
    inner_sub = AssignSub(
        Address(x, Address(Li, p)),
        Multiply(
            Address(Lx, p),
            Address(x, j)
        )
    )
    inner_loop = ForLoop(
        [
            Assign(p, Plus(Address(Lp, j), Int(1))),
            Less(p, Address(Lp, Plus(j, Int(1)))),
            PostfixInc(p)
        ],
        [Statement(inner_sub)]
    )
    outer_divide = AssignDivide(
        Address(x, j),
        Address(Lx, Address(Lp, j))
    )
    outer_loop = ForLoop(
        [
            Assign(j, Int(0)),
            Less(j, n),
            PostfixInc(j)
        ],
        [
            Statement(outer_divide),
            inner_loop
        ]
    )
    return outer_loop, inner_loop


def reachFromFiles(apath, bpath):
    A = mmread(apath)
    b = mmread(bpath)
    return reachset(A, b)