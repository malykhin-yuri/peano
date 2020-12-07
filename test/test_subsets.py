import unittest

from peano.fast_fractions import FastFraction as FF
from peano.subsets import Point, HyperFaceDivSubset
from peano.base_maps import BaseMap


class TestSubset(unittest.TestCase):

    def test_intcubes(self):
        pt = Point((FF(0, 1), FF(1, 2), FF(1, 3)))
        assert len(list(pt.gen_integer_cubes())) == 2

    def test_facediv(self):
        oc1 = HyperFaceDivSubset(dim=2, div=3, face=(0, 0), cubes=[(2,)])  # x0=0, x1>2/3
        oc2 = HyperFaceDivSubset(dim=2, div=3, face=(0, 1), cubes=[(2,)])
        bms = [bm for bm in BaseMap.gen_base_maps(dim=2, time_rev=False) if bm * oc1 == oc2]
        assert len(bms) == 1

        oc3 = HyperFaceDivSubset(dim=2, div=3, face=(0, 0))  # x0=0
        oc4 = HyperFaceDivSubset(dim=2, div=3, face=(0, 1))
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
