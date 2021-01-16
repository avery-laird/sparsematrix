from opt_ast import *
import networkx as nx


class Traversal:
    pass


class Prune(Traversal):
    """
    given a top-level for-loop and reach set,
    return two statements:
        1 the allocated reachset
        2 the pruned for-loop
    """

    def __init__(self, loop: ForLoop, reach_set: nx.DiGraph):
        self.loop = loop
        self.reach_set = reach_set

    def run(self):
        px = Iden("px")
        reachSetSize = Int(self.reach_set.order())
        reachSet = Iden("reachSet")
        # allocate space for reachSet and initialize it
        init = Statement(
            ArrayInit(IntType(), reachSet, [str(n-1) for n in nx.topological_sort(self.reach_set)])
        )
        # prune space to reachSet
        self.loop.condition[0].lhs = px
        self.loop.condition[1].lhs = px
        self.loop.condition[1].rhs = reachSetSize
        self.loop.condition[2].lhs = px
        # set value of j
        self.loop.statements.appendleft(
            Statement(Assign(Iden("j"), Address(reachSet, px)))
        )

        return init, self.loop
