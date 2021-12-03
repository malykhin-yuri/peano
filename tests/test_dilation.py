import unittest
import logging

from quicktions import Fraction

from peano import utils
from peano.dilation import Estimator

from .examples import *


def _check_dilation(data, fuzzy=False):
    curve = data['curve']
    dilation_dict = data['dilation']
    for metric in sorted(dilation_dict.keys()):
        dilation = dilation_dict[metric]
        if isinstance(dilation, list) or isinstance(dilation, tuple):
            dilation_lo, dilation_up = dilation
            dilation_eq = None
        else:
            dilation_lo = dilation * 0.999
            dilation_up = dilation * 1.001
            dilation_eq = Fraction(dilation)

        if metric == 'l2':
            func = utils.ratio_l2_squared
            dilation_lo, dilation_up = dilation_lo ** 2, dilation_up ** 2
            if dilation_eq is not None:
                dilation_eq = dilation_eq**2
        elif metric == 'l1':
            func = utils.ratio_l1
        elif metric == 'linf':
            func = utils.ratio_linf
        elif metric == 'l2_squared':
            func = utils.ratio_l2_squared

        if fuzzy:
            res = Estimator(func).estimate_dilation_fuzzy(curve, rel_tol_inv=10 ** 5)
        else:
            res = Estimator(func).estimate_dilation_regular(curve, rel_tol_inv=10 ** 5)
            if dilation_eq is not None:
                res = Estimator(func).estimate_dilation_regular(
                    curve, use_face_moments=True, face_dim=data.get('face_dim', 0), max_depth=5, rel_tol_inv=10 ** 5,
                )
                assert res['lo'] == dilation_eq
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
                'dilation': {'l2': 8, 'l1': Fraction(32, 3), 'linf': 8},
            },
            {
                'curve': get_tokarev_curve(),
                'dilation': {'l1': [98.2, 98.4], 'l2': [26.1, 26.3], 'linf': [24.1, 24.3]},
            },
            {
                'curve': get_scepin_bauman_curve(),
                'dilation': {'l1': Fraction(32, 3), 'l2': Fraction(17, 3), 'linf': Fraction(16, 3)},
            },
            {
                'curve': get_meurthe_curve(),
                'dilation': {'l1': Fraction(32, 3), 'l2': Fraction(17, 3), 'linf': Fraction(16, 3)},
            },
            {
                'curve': get_serpentine_curve(),
                'dilation': {'l1': 10, 'l2': Fraction(25, 4), 'linf': Fraction(45, 8)},
            },
            {
                'curve': get_coil_curve(),
                'dilation': {'l1': Fraction(32, 3), 'l2': Fraction(20, 3), 'linf': Fraction(20, 3)},
            },
            {
                'curve': get_r_curve(),
                'dilation': {'l1': Fraction(32, 3), 'l2': Fraction(20, 3), 'linf': Fraction(20, 3)},
            },
            {   
                'curve': get_haverkort_curve_a26(),
                'dilation': {'l1': Fraction(99 * 9 + 5, 9), 'l2': [22.7,22.9], 'linf': Fraction(12*9 + 4, 9)},
            },
            {   
                'curve': get_haverkort_curve_f(),
                'dilation': {'l1': [89.7, 89.8], 'l2': [18,19], 'linf': 14},
            },
            {
                'curve': get_ye_curve(),
                'dilation': {'l2': Fraction(408, 73)},
            },
            {
                'curve': get_spring_curve(),
                'dilation': {'l2': [16.9, 17.0]},
            },
            # Scepin & Korneev
            {
                'curve': get_tokarev_curve(),
                'dilation': {'linf': Fraction(896, 37)},
                'face_dim': 2,
            },
            {
                'curve': get_tokarev_curve(),
                'dilation': {'l2_squared': Fraction(5215408884, 7579009)},
                'face_dim': 0,
            },
        ]
        for data in known_bounds:
            _check_dilation(data)

    def test_polycurve_dilation(self):
        known_bounds = [
            # 3D curves from Haverkort Inventory
            {
                'curve': get_beta_omega_curve(),
                'dilation': {'l2': 5, 'l1': 9, 'linf': 5},
            },
            {
                'curve': get_neptunus_curve(),
                'dilation': {'l1': [88.8, 89.0], 'l2': [18.2, 18.4], 'linf': Fraction(945, 100)},
            },
            {
                'curve': get_luna_curve(),
                'dilation': {'l1': [75.5, 75.7], 'l2': [18.2, 18.4], 'linf': 14},
            },
            {
                'curve': get_iupiter_curve(),
                'dilation': {'l1': [88.6, 88.8], 'l2': [24.8, 30.0], 'linf': [16.9, 17.1]},
            },
            {
                'curve': get_ARW_Curve(),
                'dilation': {'l1': 12, 'l2': Fraction(260, 43), 'linf': Fraction(27, 5)},  # l2,l2inf: Haverkort & Walderveen
            },
        ]
        for data in known_bounds:
            _check_dilation(data)


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
