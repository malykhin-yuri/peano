import unittest

from peano.sat_adapters import CurveSATAdapter

from .examples import *  # TODO get rid of *


def get_model_from_curve(adapter, curve):
    for pnum, cnum, sp in curve.gen_defined_specs():
        adapter.add_spec_clause(pnum=pnum, cnum=cnum, spec=sp)

    if not adapter.solve():
        raise Exception("Can't get model, no such curve!")
    return adapter.get_model()


class TestSAT(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve().forget(allow_time_rev=False),
            get_meurthe_curve().forget(allow_time_rev=False),
        ]

    def test_sat(self):
        for pcurve in self.curves:
            for curve in pcurve.gen_possible_curves():
                adapter = CurveSATAdapter(pcurve)
                model = get_model_from_curve(adapter, curve)
                juncs = list(curve.gen_regular_junctions())
                for junc in juncs:
                    junc_var = adapter._get_junc_var(junc)
                    if not model[junc_var]:
                        raise Exception("Bad junc_var: False for existent junc!")
                print('.', end='', flush=True)

            print('*', flush=True)
