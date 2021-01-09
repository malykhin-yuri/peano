import unittest

from peano.utils import combinations_product


class TestUtils(unittest.TestCase):

    def test_combprod(self):
        iters = {
            1: 'abcd',
            2: 'UVW'
        }
        print(list(combinations_product([1,2,1], iters)))
