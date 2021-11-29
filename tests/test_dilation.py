import unittest

from quicktions import Fraction

from peano import utils
from peano.dilation import Estimator

from .examples import *


def _check_dilation(data, strict=False, fuzzy=False):
    curve = data['curve']
    dilation_dict = data['dilation']
    for metric in sorted(dilation_dict.keys()):
        dilation = dilation_dict[metric]
        if isinstance(dilation, list) or isinstance(dilation, tuple):
            dilation_lo, dilation_up = dilation
        else:
            dilation_lo = dilation * 0.999
            dilation_up = dilation * 1.001
            if strict:
                dilation_lo = dilation

        if metric == 'l2':
            func = utils.ratio_l2_squared
            dilation_lo, dilation_up = dilation_lo ** 2, dilation_up ** 2
        elif metric == 'l1':
            func = utils.ratio_l1
        elif metric == 'linf':
            func = utils.ratio_linf

        if fuzzy:
            res = Estimator(func).estimate_dilation_fuzzy(curve, rel_tol_inv=10 ** 5)
        elif strict:
            res = Estimator(func).estimate_dilation_regular(
                curve, use_vertex_moments=True, max_depth=data['max_depth'], rel_tol_inv=10 ** 5,
            )
            assert res['lo'] == dilation_lo
        else:
            res = Estimator(func).estimate_dilation_regular(curve, rel_tol_inv=10 ** 5)
        print(res)
        assert float(res['up']) <= dilation_up, 'metric {} up failed: {} > {}'.format(metric, res['up'], dilation_up)
        assert float(res['lo']) >= dilation_lo, 'metric {} lo failed: {} < {}'.format(metric, res['lo'], dilation_lo)


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
                'curve': get_spring_curve(),
                'dilation': {'l2': [16.9, 17.0]},
            },
        ]
        for data in known_bounds:
            _check_dilation(data)

    def test_polycurve_dilation(self):
        known_bounds = [
            {
                'curve': get_beta_omega_curve(),
                'dilation': {'l2': 5, 'l1': 9, 'linf': 5},
            },
            {
                'curve': get_neptunus_curve(),
                'dilation': {'l2': [18.2, 18.4], 'linf': [9.44, 9.46]},
            },
        ]
        for data in known_bounds:
            _check_dilation(data)

    def test_dilation_eq(self):
        known = [
            {
                'curve': get_hilbert_curve(),
                'max_depth': 10,  # TODO move to metric
                'dilation': {
                    'l2': 6,
                },
            },
            {
                'curve': get_ye_curve(),
                'max_depth': 10,
                'dilation': {'l2': Fraction(408, 73)},
            },
        ]
        for data in known:
            _check_dilation(data, strict=True)


    def test_fuzzy_dilation(self):
        known_bounds = [
            {
                'curve': get_ye_curve().forget(),
                'dilation': {'l2': [5.588, 5.59]},
            },
            {
                'curve': get_spring_curve().forget(),
                'dilation': {'l2': [16.9, 17.0]},
            },
        ]

        # special example that showed a bug
        special_curve = Curve.parse([
            ('jiJkjIJ', '0JKI~,0IKJ~,0iKJ,0Kij~,0kij,1ikJ~,0IkJ,0JIk'),
            ('jiJkjIJ', '0JKI~,0IKJ~,1iKJ,0Kij~,1kij,0ikJ~,0IkJ,0JIk'),
        ])
        known_bounds.append({
            'curve': special_curve.forget(),
            'dilation': {'l2': [18, 19]}
        })
        for data in known_bounds:
            _check_dilation(data, fuzzy=True)
