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
            ArrayInit(IntType(), reachSet, [Int(n - 1) for n in nx.topological_sort(self.reach_set)])
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


class ParallelizeOutput:
    def __init__(self,
                 level_set_init,
                 level_set_length,
                 level_set_element_length_init,
                 for_loop):
        self.level_set_init = level_set_init
        self.level_set_length = level_set_length
        self.level_set_element_length_init = level_set_element_length_init
        self.for_loop = for_loop


class Parallelize(Traversal):
    def __init__(self, loop: ForLoop, level_order_sets):
        self.loop = loop
        self.level_order_sets = level_order_sets

    def run(self) -> ParallelizeOutput:
        """
        modify the solver to execute in
        parallel, using the level order
        set to determine which loops can
        execute independently
        :param reachset:
        :return:
        """
        level_set_init = Statement(ContainerInit(VectorType(VectorType(IntType())), Iden("level_sets"),
                                             [Block(CommaSeparated([Int(e) for e in elem])) for elem in
                                              self.level_order_sets]))
        level_set_length = Statement(Assign(Iden("level_set_len"), Int(len(self.level_order_sets))))
        level_set_element_length_init = Statement(ArrayInit(IntType(), Iden("level_set_element_length"),
                                                            [Int(len(e)) for e in self.level_order_sets]))

        self.loop.condition[0].lhs = Iden("px")
        self.loop.condition[2].lhs = Iden("px")
        self.loop.condition[1] = Less(
            Iden("px"),
            Address(Iden("level_set_element_length"), Iden("i"))
        )
        self.loop.statements.appendleft(
            Statement(Assign(Iden("j"),
                             Address(Iden("level_sets"), [Iden("i"), Iden("px")])))
        )
        inner = self.loop.statements[-1]
        inner.statements.appendleft(
            Pragma("omp atomic")
        )
        return ParallelizeOutput(
            level_set_init,
            level_set_length,
            level_set_element_length_init,
            self.loop
        )
