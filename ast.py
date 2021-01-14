from enum import Enum

class Node:
    pass

class Declaration(Node):
    pass

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

class AssignDivide(Node):
    def __init__(self):
        pass

class Assign(Node):
    def __init__(self):
        pass

class AssignSub(Node):
    def __init__(self):
        pass

class Multiply(Node):
    def __init__(self, left, right):
        pass

class Increment(Node):
    def __init__(self):
        pass

class ReturnType(E)