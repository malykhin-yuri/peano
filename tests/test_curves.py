import unittest

from quicktions import Fraction

from peano.base_maps import BaseMap, Spec
from peano.curves import Curve, PathFuzzyCurve
from peano.subsets import Link
from peano.paths import PathsGenerator

from .examples import *


# some additional curves for testing
def get_rev_curves():
    chain = 'jiJ'
    bases_list = [
        'ji,Ij~,ij,JI',  # time rev at the middle
        'Ji~,ij,ij,JI',  # time rev at first cube
        'ji,ij,ij,jI~',  # time rev at last cube
    ]
    return [Curve.parse([(chain, b)]) for b in bases_list]


class TestCommon(unittest.TestCase):
    """Tests for any curves (fuzzy, poly, 2D, 3D, ...)."""

    def setUp(self):
        self.curves = [
            get_hilbert_curve(),
            get_peano_curve(),
            get_peano5_curve(),
            get_tokarev_curve(),
            get_meurthe_curve(),
            get_coil_curve(),
            get_serpentine_curve(),
            get_r_curve(),
            get_haverkort_curve_a26(),
            get_haverkort_curve_f(),

            get_hilbert_curve().forget(allow_time_rev=False),
            get_haverkort_curve_a26().forget(allow_time_rev=False),

            get_beta_omega_curve(),
            get_ARW_Curve(),
            get_neptunus_curve(),
            get_luna_curve(),
        ]
        self.curves += get_rev_curves()

    def test_rrev(self):
        """Double reverse does not change curve."""
        for curve in self.curves:
            rrcurve = ~(~curve)
            self.assertEqual(curve.proto, rrcurve.proto)
            self.assertEqual(curve.specs, rrcurve.specs)

    def test_apply_base_map_unit(self):
        """Test that the application of a sequence of base_maps with unit product does not change the curve."""
        rot_90 = BaseMap([(1,True), (0,False)])  # поворот на 90 градусов
        rot_90_t = BaseMap.parse('jI~')
        bms_list = [
            [rot_90] * 3,
            [rot_90_t] * 3,
            [
                BaseMap([(1,False), (0,False)], time_rev=True),
                BaseMap([(0,True), (1,False)]),
            ],
            [
                BaseMap([(0,True),(2,False),(1,True)], time_rev=True),
                BaseMap([(1,False),(0,True),(2,True)], time_rev=True),
            ],
        ]
        for bms in bms_list:
            # will make bms[0] * bms[1] * ... * bms[-1] * last_map = id <=> last_map = bms[-1]^{-1} * ... * bms[0]^{-1}
            last_map = BaseMap.id_map(dim=bms[0].dim)
            for bm in bms:
                last_map = bm**(-1) * last_map

            for curve in self.curves:
                if curve.dim != bms[0].dim:
                    continue
                orig = curve
                current = curve
                for bm in reversed(bms + [last_map]):
                    current = bm * current
                self.assertEqual(orig.proto, current.proto)
                self.assertEqual(orig.specs, current.specs)


class TestCurve(unittest.TestCase):
    """Tests for regular curves (peano.Curve)."""

    def setUp(self):
        self.curves = [
            get_hilbert_curve(),
            get_peano_curve(),
            get_peano5_curve(),
            get_scepin_bauman_curve(),
            get_tokarev_curve(),
            get_meurthe_curve(),
            get_coil_curve(),
            get_serpentine_curve(),
            get_r_curve(),
            get_haverkort_curve_a26(),
            get_haverkort_curve_f(),
            get_beta_omega_curve(),
            get_ARW_Curve(),
            get_neptunus_curve(),
            get_luna_curve(),
            get_spring_curve(),
        ]
        self.curves += get_rev_curves()

    def test_check(self):
        for curve in self.curves:
            curve.check()
            (~curve).check()

    def test_fractions(self):
        for curve in self.curves:
            for cnum in range(curve.genus):
                fraction = curve.specs[cnum] * curve
                fraction.check()

    def test_subdivision(self):
        for curve in self.curves:
            curve.get_subdivision().check()

    def test_subsubdivision(self):
        for curve in self.curves:
            curve.get_subdivision(2).check()

    def test_compose(self):
        for curve in self.curves:
            for bm in BaseMap.gen_base_maps(dim=curve.dim):
                for pnum in range(curve.pcount):
                    for cnum in range(curve.genus):
                        spec = Spec(base_map=bm, pnum=pnum)
                        C = spec * curve
                        assert C.specs[cnum] * C == curve._compose_specs(spec, cnum) * curve

    def test_junc(self):
        for curve in self.curves:
            for jcnt, junc in enumerate(curve.gen_regular_junctions()):
                if jcnt > 100:
                    raise Exception("Too many juncs!")
            self.assertEqual(set(curve.gen_regular_junctions()), set(curve.get_junctions_info().keys()))

    def test_vertex_moments(self):
        known_moments = [
            {
                'curve': get_haverkort_curve_a26(),
                'moments': [Fraction(k, 28) for k in [0, 5, 9, 12, 16, 19, 23, 28]],
            },
            {
                'curve': get_haverkort_curve_f(),
                'moments': [Fraction(k, 28) for k in [1, 6, 8, 13, 15, 20, 22, 27]],
            },
            {
                'curve': get_tokarev_curve(),
                'moments': [Fraction(k, 126) for k in [0, 22, 41, 50, 76, 85, 104, 126]],
            },
        ]
        for d in known_moments:
            moments = d['curve'].get_vertex_moments().values()
            self.assertSetEqual(set(d['moments']), set(moments))

    def test_gate(self):
        known = [
            {
                'curve': get_haverkort_curve_f(),
                'gates': [Link.parse_gates('(0,1/3,1/3)->(2/3,1/3,0)')],
            },
            {
                'curve': get_beta_omega_curve(),
                'gates': [Link.parse_gates('(0,1/3)->(1,1/3)'), Link.parse_gates('(0,1/3)->(2/3,0)')],
            },
            {
                'curve': get_neptunus_curve(),
                'gates': [Link.parse_gates('(0,0,0)->(1,0,0)'), Link.parse_gates('(0,0,0)->(1,1,1)')],
            },
            {
                'curve': get_luna_curve(),
                'gates': [Link.parse_gates('(0,0,0)->(1,0,0)'), Link.parse_gates('(0,0,0)->(1,1,1)')],
            },
        ]
        for data in known:
            curve = data['curve']
            gates = [Link(curve.get_entrance(pnum), curve.get_exit(pnum)) for pnum in range(curve.pcount)]
            self.assertEqual(gates, data['gates'])

            for pnum, gate in enumerate(gates):
                entr_face = [1 if pj == Fraction(1) else 0 if pj == Fraction(0) else None for pj in gate.entrance]
                assert curve.get_face_moment(entr_face, pnum) == Fraction(0)

                exit_face = [1 if pj == Fraction(1) else 0 if pj == Fraction(0) else None for pj in gate.exit]
                assert curve.get_face_moment(exit_face, pnum, last=True) == Fraction(1)

    def test_depth(self):
        known = [
            {'curve': get_hilbert_curve(), 'depth': 2},
            {'curve': get_peano_curve(), 'depth': 1},
            {'curve': get_tokarev_curve(), 'depth': 3},
        ]
        for data in known:
            assert data['curve'].get_depth() == data['depth']


class TestFuzzyCurves(unittest.TestCase):
    """Tests specific for fuzzy curves."""

    def setUp(self):
        self.curves = [
            get_hilbert_curve(),
            get_peano_curve(),
            get_tokarev_curve(),
        ]

    def test_check(self):
        for curve in [get_hilbert_curve()]:
            pcurve = curve.forget(allow_time_rev=True)
            for c in pcurve.gen_possible_curves():
                c.check()

    def test_junc(self):

        def is_specialization(curve, tmpl):
            # Check if curve has all of defined specs of the given template, and they are the same."""
            return all(curve.patterns[pnum].specs[cnum] == sp for pnum, cnum, sp in tmpl.gen_defined_specs())

        for num, curve in enumerate(self.curves):
            if num == 0:
                pcurve = curve.forget(allow_time_rev=True)
            else:
                pcurve = curve.forget(allow_time_rev=False)

            junc_info = pcurve.get_junctions_info()
            for junc in junc_info:
                if any(dx not in set([0,1,-1]) for dx in junc.delta_x):
                    raise Exception("Bad delta!")

            for curve in pcurve.gen_possible_curves():
                juncs = list(curve.gen_regular_junctions())
                # проверяем, что для каждого найденного стыка есть порождающая кривая
                for junc in juncs:
                    if junc not in junc_info:
                        raise Exception("Unknown junc!")
                    if not any(is_specialization(curve, tmpl) for tmpl in junc_info[junc]):
                        raise Exception("Can't found consistent curve!")

                # проверяем, что для каждого не найденного стыка нет порождающих кривых
                for junc, curves in junc_info.items():
                    if junc in juncs:
                        continue
                    if any(is_specialization(curve, tmpl) for tmpl in junc_info[junc]):
                        raise Exception("Found consistent curve for wrong junc!")


class TestMisc(unittest.TestCase):
    def test_diag(self):
        gates = Link.parse_gates('(0,0)->(1,1/2)')
        pgen = PathsGenerator(dim=2, div=5, links=[gates], max_cdist=1)
        paths = next(pgen.generate_paths())
        pcurve = PathFuzzyCurve.init_from_paths(paths)
        curve = pcurve.get_curve_example()
        print(curve)
