from enum import Enum
from collections import deque


class Node:
    def gen(self) -> str:
        raise Exception('node is abstract!')


class Type(Node):
    pass


class Statement(Node):
    def __init__(self, unit):
        self.unit = unit

    def gen(self) -> str:
        return "{};".format(self.unit.gen())


class LoopConstruct(Node):
    def __init__(self, condition: list[Node], statements: list[Statement]):
        self.condition = condition
        self.statements = deque(statements)
        self.name = None

    def gen_condition(self) -> str:
        return self.condition[0].gen()

    def gen(self) -> str:
        return "{} ({}) {{\n{}\n}}".format(self.name.gen(), self.gen_condition(),
                                           "\n".join([s.gen() for s in self.statements]))


class ForLoop(LoopConstruct):
    def __init__(self, condition: list[Node], statements: list[Statement]):
        if len(condition) > 3:
            raise Exception("at most 3 conditions allowed")
        super().__init__(condition, statements)
        self.name = Reserved("for")

    def gen_condition(self) -> str:
        return "; ".join([c.gen() for c in self.condition])

    def replace_condition(self, new_cond: Node):
        self.condition[1] = new_cond

    def prepend_statement(self, new_statment: Statement):
        self.statements.appendleft(new_statment)


class Reserved(Node):
    def __init__(self, name: str):
        self.name = name

    def gen(self) -> str:
        return self.name


class IntType(Type):
    def gen(self):
        return 'int'


class VectorType(Type):
    def __init__(self, elem_type):
        self.elem_type = elem_type

    def gen(self) -> str:
        return "std::vector<{}>".format(self.elem_type.gen())


class Number(Node):
    def __init__(self, value):
        self.value = value

    def gen(self) -> str:
        return str(self.value)


class Int(Number):
    def __init__(self, value: int):
        super().__init__(value)


class Float(Number):
    def __init__(self, value: float):
        super().__init__(value)


class Iden(Node):
    def __init__(self, name: str):
        self.name = name

    def gen(self):
        return str(self.name)


class Declaration(Node):
    def __init__(self, t: Type, val: Node = None):
        self.type = t
        self.val = val

    def gen(self):
        return "{}{};".format(self.type.gen(), "" if not self.val else " " + self.val.gen())


class ArrayInit(Node):
    """
    special class for array init, because it's
    more convenient to implement just a subset
    of a full C++ AST
    """

    def __init__(self, t: Type, name: Iden, values: list[Node]):
        self.type = t
        self.values = values
        self.name = name

    def gen(self) -> str:
        return "{} {}[] = {{{}}}".format(self.type.gen(), self.name.gen(), ",".join([n.gen() for n in self.values]))


class ContainerInit(ArrayInit):
    def gen(self) -> str:
        return "{} {} = {{{}}}".format(self.type.gen(), self.name.gen(), ",".join([n.gen() for n in self.values]))


class CommaSeparated(Node):
    def __init__(self, values: list[Node]):
        self.values = values

    def gen(self) -> str:
        return ",".join([n.gen() for n in self.values])


class Block(Node):
    """
    Represents a scoped block.
    """

    def __init__(self, value: Node):
        self.value = value

    def gen(self) -> str:
        return "{{{}}}".format(self.value.gen())


class Symbol(Node):
    pass


class Branch(Node):
    def __init__(self, true, false):
        pass


class Return(Node):
    def __init__(self):
        pass


class For(Node):
    def __init__(self):
        pass


class BinaryOp(Node):
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs
        self.symbol = None

    def gen(self) -> str:
        return "{} {} {}".format(self.lhs.gen(), self.symbol, self.rhs.gen())


class Assign(BinaryOp):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)
        self.symbol = "="


class Less(BinaryOp):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)
        self.symbol = "<"


class AssignDivide(BinaryOp):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)
        self.symbol = "/="


class AssignSub(BinaryOp):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)
        self.symbol = "-="


class Multiply(BinaryOp):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)
        self.symbol = "*"


class Plus(BinaryOp):
    def __init__(self, lhs, rhs):
        super().__init__(lhs, rhs)
        self.symbol = "+"


class PostfixOp(Node):
    def __init__(self, lhs):
        self.lhs = lhs
        self.symbol = None

    def gen(self) -> str:
        return "{}{}".format(self.lhs.gen(), self.symbol.gen())


class PostfixInc(PostfixOp):
    def __init__(self, lhs):
        super().__init__(lhs)
        self.symbol = Reserved("++")


class Address(Node):
    def __init__(self, target: Node, addr):
        self.target = target
        self.addr = addr

    def gen(self) -> str:
        if isinstance(self.addr, list):
            return "{}{}".format(self.target.gen(), "".join(["[{}]".format(x.gen()) for x in self.addr]))
        else:
            return "{}[{}]".format(self.target.gen(), self.addr.gen())


class Pragma(Node):
    """
    A very low level node for inserting preprocessor
    statements into the AST. It only stores the text
    to be printed, and nothing else about the internal
    structure.
    """

    def __init__(self, text: str):
        self.text = text

    def gen(self) -> str:
        return "#pragma {}".format(self.text)
