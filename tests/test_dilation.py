import unittest

from sympy import Rational

from peano import utils
from peano.paths import PathsGenerator
from peano.curves import PathFuzzyCurve, Proto
from peano.dilation import Estimator
from peano.base_maps import BaseMap

from .examples import *  # TODO get rid of "*"


class TestCurve(unittest.TestCase):

    def test_curve_dilation(self):
        known_bounds = [
            {
                'curve': get_hilbert_curve(),
                'dilation': { 'l2': 6, 'l1': 9, 'linf': 6},
            },
            {
                'curve': get_peano_curve(),
                'dilation': {'l2': 8, 'l1': 32/3, 'linf': 8},
            },
            {
                'curve': get_tokarev_curve(),
                'dilation': {'l1': [98.2, 98.4], 'l2': [26.1, 26.3], 'linf': [24.1, 24.3]},
            },
            {
                'curve': get_scepin_bauman_curve(),
                'dilation': {'l1': (10 + 2/3), 'l2': (5 + 2/3), 'linf': (5 + 1/3)},
            },
            {
                'curve': get_meurthe_curve(),
                'dilation': {'l1': (10 + 2/3), 'l2': (5 + 2/3), 'linf': (5 + 1/3)},
            },
            {
                'curve': get_serpentine_curve(),
                'dilation': {'l1': 10, 'l2': 6.25, 'linf': 5.625},
            },
            {
                'curve': get_coil_curve(),
                'dilation': {'l1': (10 + 2/3), 'l2': (6 + 2/3), 'linf': (6 + 2/3)},
            },
            {
                'curve': get_r_curve(),
                'dilation': {'l1': (10 + 2/3), 'l2': (6 + 2/3), 'linf': (6 + 2/3)},
            },
            {   
                'curve': get_haverkort_curve_a26(),
                'dilation': {'l1': (99 + 5/9), 'l2': [22.7,22.9], 'linf': (12 + 4/9)},
            },
            {   
                'curve': get_haverkort_curve_f(),
                'dilation': {'l1': [89.7, 89.8], 'l2': [18,19], 'linf': 14},
            },
            {
                'curve': get_ye_curve(),
                'dilation': {'l2': [5.588, 5.59]},
            },
            {
                'curve': get_17_curve(),
                'dilation': {'l2': [16.9, 17.0]},
            },
        ]
        for data in known_bounds:
            curve = data['curve']
            dilation_dict = data['dilation']
            for metric in sorted(dilation_dict.keys()):
                dilation = dilation_dict[metric]
                if isinstance(dilation, list):
                    dilation_lo, dilation_up = dilation
                else:
                    dilation_lo = dilation * 0.999
                    dilation_up = dilation * 1.001

                if metric == 'l2':
                    func = utils.ratio_l2_squared
                    dilation_lo, dilation_up = dilation_lo**2, dilation_up**2
                elif metric == 'l1':
                    func = utils.ratio_l1
                elif metric == 'linf':
                    func = utils.ratio_linf

                res = Estimator(func).estimate_dilation(curve, rel_tol_inv=10 ** 5, verbose=False)
                print(res)
                assert float(res['up']) <= dilation_up, 'metric {} up failed: {} > {}'.format(metric, res['up'], dilation_up)
                assert float(res['lo']) >= dilation_lo, 'metric {} lo failed: {} < {}'.format(metric, res['lo'], dilation_lo)

    def test_polycurve_dilation(self):
        known_bounds = [
            {
                'curve': get_beta_omega_curve(),
                'dilation': { 'l2': 5, 'l1': 9, 'linf': 5},
            },
            {
                'curve': get_neptunus_curve(),
                'dilation': { 'l2': [18.2, 18.4], 'linf': [9.44, 9.46]},
            },
        ]
        for data in known_bounds:
            curve = data['curve']
            for metric, dilation in data['dilation'].items():
                if isinstance(dilation, list):
                    dilation_lo, dilation_up = dilation
                else:
                    dilation_lo = dilation * 0.999
                    dilation_up = dilation * 1.001

                if metric == 'l2':
                    func = utils.ratio_l2_squared
                    dilation_lo, dilation_up = dilation_lo**2, dilation_up**2
                elif metric == 'l1':
                    func = utils.ratio_l1
                elif metric == 'linf':
                    func = utils.ratio_linf

                res = Estimator(func).estimate_dilation(curve, rel_tol_inv=10 ** 4)
                assert float(res['up']) <= dilation_up, 'metric {} up failed: {} > {}'.format(metric, res['up'], dilation_up)
                assert float(res['lo']) >= dilation_lo, 'metric {} lo failed: {} < {}'.format(metric, res['lo'], dilation_lo)

    def test_ye_dilation(self):
        # TODO: use ye curve from examples ?
        good_proto = Proto(dim=2, div=5, cubes=[
            (0, 0), (0, 1), (1, 1), (1, 0), (2, 0),
            (2, 1), (2, 2), (1, 2), (0, 2), (0, 3),
            (0, 4), (1, 4), (1, 3), (2, 3), (2, 4),
            (3, 4), (4, 4), (4, 3), (3, 3), (3, 2),
            (4, 2), (4, 1), (3, 1), (3, 0), (4, 0),
        ])
        # in new version we have (0,0)->(0,1) gate
        good_proto = BaseMap.parse('ji') * good_proto

        paths_gen = PathsGenerator(dim=2, div=5, hdist=1, max_cdist=1)
        for paths in paths_gen.generate_paths():
            if paths[0].proto == good_proto:
                path0 = paths[0]
                break

        pcurve = PathFuzzyCurve.init_from_paths([path0])
        estimator = Estimator(utils.ratio_l2_squared)
        curve = estimator.estimate_dilation(pcurve, rel_tol_inv=10000)['curve']
        dilation = estimator.estimate_dilation(curve, rel_tol_inv=None, use_vertex_brkline=True, max_depth=5)

        assert dilation['lo'] == (Rational(408, 73) ** 2)
