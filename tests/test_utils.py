import unittest

from peano.utils import combinations_product


class TestUtils(unittest.TestCase):

    def test_combprod(self):
        iters = {
            1: 'abcd',
            2: 'UVW'
        }
        output = list(combinations_product([1,2,1], iters))
        print(output)
        self.assertEqual(len(output), 30)
        self.assertTrue(('a','V','d') in output)
        self.assertTrue(('b','U','a') not in output)
