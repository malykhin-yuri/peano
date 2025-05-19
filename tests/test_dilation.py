import unittest
import logging
from typing import Any
from dataclasses import asdict

from quicktions import Fraction  # type: ignore

from peano import utils
from peano.dilation import Estimator
from peano.zoo import CurveInfo, DilationInfo, get_all_curves, get_ye_curve, get_spring_curve
from peano.curves import Curve


def _check_dilation(data: CurveInfo, fuzzy: bool = False) -> None:
    curve = data.curve
    dilation_dict = asdict(data.dilation)
    for metric in sorted(dilation_dict.keys()):
        if metric == 'spec':
            continue  # UGLY TODO FIX
        dilation = dilation_dict[metric]
        if dilation is None:
            continue
        if isinstance(dilation, list) or isinstance(dilation, tuple):
            dilation_lo, dilation_up = dilation
            dilation_eq = None
        else:
            if isinstance(dilation, dict):
                print('ZOPA BLIAT', data.dilation, metric)
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
                if data.dilation.spec is not None:
                    face_dim = data.dilation.spec.get(metric, {}).get('face_dim', 0)
                else:
                    face_dim = 0
                res = Estimator(func).estimate_dilation_regular(
                    curve, use_face_moments=True, face_dim=face_dim, max_depth=5, rel_tol_inv=10 ** 5,  # WHY max_depth = 5?
                )
                assert res['lo'] == dilation_eq
        print(res)
        assert float(res['up']) <= dilation_up, 'metric {} up failed: {} > {}'.format(metric, res['up'], dilation_up)
        assert float(res['lo']) >= dilation_lo, 'metric {} lo failed: {} < {}'.format(metric, res['lo'], dilation_lo)


class TestCurve(unittest.TestCase):

    def test_curve_dilation(self) -> None:
        for data in get_all_curves():
            if data.dilation is not None:
                _check_dilation(data)


    def test_fuzzy_dilation(self) -> None:
        known_bounds = [
            CurveInfo(
                curve=get_ye_curve().curve.forget(),
                dilation=DilationInfo(l2=[5.588, 5.59]),
            ),
            CurveInfo(
                curve=get_spring_curve().curve.forget(),
                dilation=DilationInfo(l2=[16.9, 17.0]),
            ),
        ]

        # special example that showed a bug
        special_curve = Curve.parse([
            ('jiJkjIJ', '0JKI~,0IKJ~,0iKJ,0Kij~,0kij,1ikJ~,0IkJ,0JIk'),
            ('jiJkjIJ', '0JKI~,0IKJ~,1iKJ,0Kij~,1kij,0ikJ~,0IkJ,0JIk'),
        ])
        known_bounds.append(CurveInfo(
            curve=special_curve.forget(),
            dilation=DilationInfo(l2=[18, 19]),
        ))
        for data in known_bounds:
            _check_dilation(data, fuzzy=True)
