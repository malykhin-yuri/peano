import unittest

from peano.sat_adapters import CurveSATAdapter

from .examples import *  # TODO get rid of *


def get_model_from_curve(adapter, curve):
    for pnum, cnum, sp in curve.sp_info():
        sp_var = adapter.get_sp_var(pnum, cnum, sp)
        adapter.append_clause({sp_var: True})

    if not adapter.solve():
        raise Exception("Can't get model, no such curve!")
    return adapter.get_model()


class TestSAT(unittest.TestCase):
    def setUp(self):
        self.curves = [
            get_peano_curve().forget(),
            get_meurthe_curve().forget(),
        ]

    def test_sat(self):
        for pcurve in self.curves:
            for curve in pcurve.gen_possible_curves():
                adapter = CurveSATAdapter(dim=pcurve.dim)
                adapter.init_curve(pcurve)
                model = get_model_from_curve(adapter, curve)
                juncs = list(curve.gen_regular_junctions())
                for junc in juncs:
                    junc_var = adapter.get_junc_var(junc)
                    if not model[junc_var]:
                        raise Exception("Bad junc_var: False for existent junc!")
                print('.', end='', flush=True)

            print('*', flush=True)
