import unittest

from quicktions import Fraction

from peano.subsets import Point, FacetDivSubset, Link
from peano.base_maps import BaseMap


class TestSubset(unittest.TestCase):

    def test_intcubes(self):
        pt = Point((Fraction(0, 1), Fraction(1, 2), Fraction(1, 3)))
        assert len(list(pt._gen_integer_cubes())) == 2

        pt = Point((Fraction(0), Fraction(0), Fraction(0)))
        assert len(list(pt._gen_integer_cubes())) == 8

    def test_facediv(self):
        oc1 = FacetDivSubset(dim=2, div=3, facet=(0, 0), cubes=[(2,)])  # x0=0, x1>2/3
        oc2 = FacetDivSubset(dim=2, div=3, facet=(0, 1), cubes=[(2,)])
        bms = [bm for bm in BaseMap.gen_base_maps(dim=2, time_rev=False) if bm * oc1 == oc2]
        assert len(bms) == 1

        oc3 = FacetDivSubset(dim=2, div=3, facet=(0, 0))  # x0=0
        oc4 = FacetDivSubset(dim=2, div=3, facet=(0, 1))
        bms = [bm for bm in BaseMap.gen_base_maps(dim=2, time_rev=False) if bm * oc3 == oc4]
        assert len(bms) == 2

        assert oc3.map_to_cube(div=3, cube=(0, 2)) == oc1
        assert oc4.map_to_cube(div=3, cube=(2, 2)) == oc2

        assert len(list(oc1.gen_neighbours())) == 1
        assert list(oc1.gen_neighbours()) == [((-1, 0), oc2)]

        assert len(list(oc1.divide(div=3))) == 1
        assert len(list(oc2.divide(div=3))) == 1

        assert len(list(oc3.divide(div=3))) == 3
        assert len(list(oc4.divide(div=3))) == 3

        oc5 = FacetDivSubset(dim=2, div=5, facet=(1, 1), cubes=[(2,)])
        oc5_div = list(oc5.divide(5))
        self.assertEqual(oc5_div[0][0], (2, 4))
        self.assertEqual(oc5_div[0][1], FacetDivSubset(dim=2, div=5, facet=(1, 1), cubes=[]))

    def test_std(self):
        gates = Link.parse_gates('(1,1/2,0)->(0,1/3,2/3)')
        print(gates.std())
