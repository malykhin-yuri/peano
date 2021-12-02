import unittest

from peano._sat_adapters import CurveSATAdapter
from peano.curves import Curve

from .examples import get_peano_curve, get_ye_curve


class TestSAT(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve(),
            get_ye_curve(),
        ]

    def test_sat(self):
        for curve in self.curves:
            adapter = CurveSATAdapter(curve=curve)
            adapter.solve()
            model_curve = adapter.get_model_curve()
            self.assertEqual(curve, model_curve)

            empty_patterns = []
            for pattern in curve.patterns:
                empty_pattern = (pattern.proto, (None,) * curve.genus)
            empty_curve = Curve(dim=curve.dim, div=curve.div, patterns=empty_patterns)

            for junc in curve.gen_regular_junctions():
                adapter2 = CurveSATAdapter(curve=curve)
                adapter2.add_forbid_clause(junc, empty_curve)
                self.assertFalse(adapter2.solve())
