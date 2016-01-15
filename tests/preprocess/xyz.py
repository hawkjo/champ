import unittest
import numpy as np
from preprocess.xyz import XYZFile


class XYZFileTests(unittest.TestCase):
    def test_str(self):
        f = XYZFile(np.array([[1, 2, 3],
                              [4, 5, 6]]))
        expected = """0 0 1
1 0 2
2 0 3
0 1 4
1 1 5
2 1 6
"""
        self.assertEqual(str(f), expected)
