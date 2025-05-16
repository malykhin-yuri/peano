import unittest

from quicktions import Fraction  # type: ignore

from peano.utils import combinations_product, get_periodic_sum


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

    def test_periodic_sum(self):
        self.assertEqual(get_periodic_sum([1, 2, 3], [4, 5], 6), Fraction(907, 3780))
