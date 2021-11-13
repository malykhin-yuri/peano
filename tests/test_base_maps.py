import unittest

from quicktions import Fraction

from peano.subsets import Point
from peano.base_maps import BaseMap


class TestBaseMap(unittest.TestCase):
    def setUp(self):
        self.base_maps = []
        for dim in range(2, 6):
            self.base_maps += list(BaseMap.gen_base_maps(dim))

    def test_mul(self):
        bm1 = BaseMap.parse('Ij')
        bm2 = BaseMap.parse('ji')
        self.assertEqual(bm1 * bm2, BaseMap.parse('jI'))
        self.assertEqual(bm2 * bm1, BaseMap.parse('Ji'))

    def test_inv(self):
        for bm in self.base_maps:
            self.assertEqual(bm * bm**(-1), BaseMap.id_map(dim=bm.dim))
            self.assertEqual(bm**(-1) * bm, BaseMap.id_map(dim=bm.dim))

    def test_conj(self):
        dim3 = [bm for bm in self.base_maps if bm.dim == 3]
        for bm1 in dim3:
            for bm2 in dim3:
                # TODO - правилен ли такой конжугад?
                self.assertEqual(bm1.conjugate_by(bm2), bm2 * bm1 * bm2**(-1))

    def test_constraint_fast(self):
        pairs = [
            (
                [Fraction(0, 1), Fraction(1, 1), Fraction(1, 3), Fraction(1, 2)],
                [Fraction(0, 1), Fraction(1, 2), Fraction(1, 3), Fraction(0, 1)],
            ),
            (
                [Fraction(0, 1), Fraction(0, 1)],
                [Fraction(1, 2), Fraction(1, 2)],
            ),
            (
                [Fraction(0, 1), Fraction(0, 1), Fraction(1, 2), Fraction(1, 2), Fraction(1, 4)],
                [Fraction(1, 2), Fraction(1, 2), Fraction(1, 1), Fraction(1, 1), Fraction(3, 4)],
            ),
            (
                [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)],
                [Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(1, 1)],
            ),
        ]
        for src_coords, dst_coords in pairs:
            src = Point(src_coords)
            dst = Point(dst_coords)
            dim = len(src)
            bms = set(BaseMap.gen_constraint_fast(src, dst))
            for bm in BaseMap.gen_base_maps(dim, time_rev=False):
                if bm in bms:
                    assert bm * src == dst
                else:
                    assert bm * src != dst
