import unittest

from peano._sat_adapters import CurveSATAdapter
from peano.curves import Curve, FuzzyCurve
from peano.paths import Proto

from peano.zoo import get_all_curves, get_hilbert_curve, get_peano_curve, get_beta_omega_curve


def get_empty_curve(curve: Curve) -> FuzzyCurve:
    empty_patterns: list[tuple[Proto, tuple[None, ...]]] = []
    for pattern in curve.patterns:
        empty_patterns.append((pattern.proto, (None,) * curve.genus))
    return FuzzyCurve(dim=curve.dim, div=curve.div, patterns=empty_patterns)


class TestSAT(unittest.TestCase):
    def setUp(self):
        self.curves = [x.curve for x in get_all_curves()]
        self.fuzzy_curves = [
            get_hilbert_curve().curve.forget(),
            get_peano_curve().curve.forget(),
            get_beta_omega_curve().curve.forget(),
        ]

    def test_has_model(self) -> None:
        for curve in self.curves:
            adapter = CurveSATAdapter(curve=curve)
            adapter.solve()
            model_curve = adapter.get_model_curve()
            self.assertEqual(curve, model_curve)

    def test_fuzzy_has_model(self) -> None:
        for fuzzy_curve in self.fuzzy_curves:
            adapter = CurveSATAdapter(curve=fuzzy_curve)
            self.assertTrue(adapter.solve())

    def test_forbidden_subcurve(self) -> None:
        # assert that if we forbid a subcurve, no model is found
        for curve in self.curves:
            juncs = list(curve.gen_auto_junctions()) + list(curve.gen_regular_junctions())
            for junc in juncs:
                for pnum in range(curve.mult):
                    for cnum in range(curve.genus):
                        subcurve = get_empty_curve(curve)
                        spec = curve.patterns[pnum].specs[cnum]
                        # TODO Can't use specify method because this is not PathFuzzyCurve :(
                        subcurve = subcurve._specify_allowed({pnum: {cnum: spec}})
                        adapter = CurveSATAdapter(curve=curve)
                        adapter.add_forbid_clause(junc, subcurve)
                        self.assertFalse(adapter.solve())

    def test_forbidden_junctions(self) -> None:
        # assert that if we forbid any existing junction, no model is found
        for curve in self.curves:
            empty_curve = get_empty_curve(curve)
            for junc in curve.gen_regular_junctions():
                adapter = CurveSATAdapter(curve=curve)
                adapter.add_forbid_clause(junc, empty_curve)
                self.assertFalse(adapter.solve())

    def test_forbidden_junctions_fuzzy(self) -> None:
        # assert that if we forbid all possible junctions, no model is found
        for fuzzy_curve in self.fuzzy_curves:
            empty_curve = get_empty_curve(fuzzy_curve)
            adapter = CurveSATAdapter(curve=fuzzy_curve)
            for junc in fuzzy_curve.gen_regular_junctions():
                adapter.add_forbid_clause(junc, empty_curve)
            self.assertFalse(adapter.solve())
