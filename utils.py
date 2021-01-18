from scipy.io import mmread
from symbolicanalysis import reachset, level_order_set
from opt_ast import *


def make_naive_solve():
    """
    This function constructs a representation of the
    naive solver included with the project specification.
    :return:
    """
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


def readMatrix(apath, bpath):
    A = mmread(apath)
    b = mmread(bpath)
    return A, b


def reachFromFiles(apath, bpath):
    A, b = readMatrix(apath, bpath)
    return reachset(A, b), A, b


def levelOrderFromFiles(apath, bpath):
    A, _ = readMatrix(apath, bpath)
    return level_order_set(A)
