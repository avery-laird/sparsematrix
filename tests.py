import unittest
import numpy as np
from scipy.sparse.linalg import spsolve_triangular
from scipy.io import mmread

class Correctness(unittest.TestCase):
    def test_answer_correct(self):
        A = mmread('torso1/torso1.mtx')
        b = mmread('b_for_torso1.mtx')
        x = spsolve_triangular(A, b)
        



if __name__ == '__main__':
    unittest.main()
