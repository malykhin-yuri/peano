import unittest
import itertools

from quicktions import Fraction  # type: ignore

from peano.base_maps import BaseMap, Spec
from peano.curves import Curve, Junction, PathFuzzyCurve
from peano.subsets import Link
from peano.paths import PathsGenerator
from peano.gates import GatesGenerator
from peano.utils import gen_faces
from peano.zoo import *


# some additional curves for testing
def get_rev_curves() -> list[Curve]:
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
        self.curves = [x.curve for x in get_all_curves()]
        self.curves.extend([
            get_hilbert_curve().curve.forget(disable_time_rev=True),
            get_haverkort_curve_a26().curve.forget(disable_time_rev=True),
        ])
        self.curves += get_rev_curves()

    def test_rrev(self) -> None:
        """Double reverse does not change curve."""
        for curve in self.curves:
            rrcurve = ~(~curve)
            self.assertEqual(curve.proto, rrcurve.proto)
            self.assertEqual(curve.specs, rrcurve.specs)

    def test_apply_base_map_unit(self) -> None:
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
        self.curves_info = get_all_curves()

        # TODO: maybe some filtering is required
        self.curves = [x.curve for x in self.curves_info]
        self.curves += get_rev_curves()

    def test_check(self) -> None:
        for curve in self.curves:
            curve.check()
            (~curve).check()

    def test_fractions(self) -> None:
        for curve in self.curves:
            for cnum in range(curve.genus):
                fraction = curve.specs[cnum] * curve
                fraction.check()

    def test_subdivision(self) -> None:
        for curve in self.curves:
            curve.get_subdivision().check()

    def test_subsubdivision(self) -> None:
        for curve in self.curves:
            curve.get_subdivision(2).check()
            self.assertEqual(curve.get_subdivision().get_subdivision(), curve.get_subdivision(3))

    def test_compose(self) -> None:
        for curve in self.curves:
            for bm in BaseMap.gen_base_maps(dim=curve.dim):
                for pnum in range(curve.mult):
                    for cnum in range(curve.genus):
                        spec = Spec(base_map=bm, pnum=pnum)
                        C = spec * curve
                        assert C.specs[cnum] * C == curve._compose_specs(spec, cnum) * curve

    def test_junc(self) -> None:
        def is_continuous(curve: Curve, junc: Junction) -> bool:
            link1 = junc.spec1.base_map * Link(curve.get_entrance(junc.spec1.pnum), curve.get_exit(junc.spec1.pnum))
            link2 = junc.spec2.base_map * Link(curve.get_entrance(junc.spec2.pnum), curve.get_exit(junc.spec2.pnum))
            return link1.exit == link2.entrance.transform(shift=junc.delta_x)

        for curve in self.curves:
            for jcnt, junc in enumerate(curve.gen_regular_junctions()):
                self.assertTrue(is_continuous(curve, junc))
                self.assertTrue(all(x in [0,-1,1] for x in junc.delta_x))
                self.assertNotEqual(junc.delta_x, (0,) * curve.dim)
                if jcnt > 100:
                    raise Exception("Too many juncs!")
            self.assertEqual(set(curve.gen_regular_junctions()), set(curve.get_junction_templates()))


    def test_gate(self):
        for info in self.curves_info:
            if info.gates is None:
                continue
            curve = info.curve
            gates = [Link(curve.get_entrance(pnum), curve.get_exit(pnum)) for pnum in range(curve.mult)]
            self.assertEqual(gates, info.gates)

            for pnum, gate in enumerate(gates):
                entr_face = [1 if pj == Fraction(1) else 0 if pj == Fraction(0) else None for pj in gate.entrance]
                assert curve.get_face_moment(entr_face, pnum) == Fraction(0)

                exit_face = [1 if pj == Fraction(1) else 0 if pj == Fraction(0) else None for pj in gate.exit]
                assert curve.get_face_moment(exit_face, pnum, last=True) == Fraction(1)

    def test_other_info(self):
        for info in self.curves_info:
            if info.moments is not None:
                for face_dim, expected_moments in info.moments.items():
                    last = (face_dim < 0)
                    got_moments = [info.curve.get_face_moment(face, last=last) for face in gen_faces(info.curve.dim, abs(face_dim))]
                    self.assertEqual(list(sorted(expected_moments)), list(sorted((got_moments))))

            if info.junctions_count is not None:
                juncs = list(info.curve.gen_regular_junctions())
                self.assertEqual(info.junctions_count, len(juncs))

            if info.depth is not None:
                self.assertEqual(info.depth, info.curve.get_depth())


class TestFuzzyCurves(unittest.TestCase):
    """Tests specific for fuzzy curves."""

    def setUp(self):
        self.curves = [
            get_hilbert_curve().curve,
            get_peano_curve().curve,
            get_tokarev_curve().curve,
        ]

    def test_check(self) -> None:
        for curve in [get_hilbert_curve().curve]:
            pcurve = curve.forget()
            for c in pcurve.gen_possible_curves():
                c.check()

    def test_GP(self) -> None:
        # Haverkort & Walderveen, p.135: "... 272 orders"""
        all_bms = list(BaseMap.gen_base_maps(dim=2))
        pg = PathsGenerator(dim=2, div=3, links=[Link.parse_gates('(0,0)->(1,1)')], max_cdist=1)
        paths_tuple = pg.get_paths_example()
        seen = set()
        fcurve = PathFuzzyCurve.init_from_paths(paths_tuple, disable_time_rev=True)
        for curve in fcurve.gen_possible_curves():
            if all(bm * curve not in seen for bm in all_bms):
                seen.add(curve)
        self.assertEqual(len(seen), 272)

    def test_3D(self) -> None:
        # Haverkort inventory: 920 face-continuous vertex-gated order-preserving
        dim, div, mult = 3, 2, 1
        all_bms = list(BaseMap.gen_base_maps(dim=dim))
        seen = set()
        gates_gen = GatesGenerator(dim=dim, div=div, mult=mult)
        for gates in gates_gen.gen_gates():
            link = gates[0]
            if link.entrance.face_dim() > 0 or link.exit.face_dim() > 0:
                continue
            paths_gen = PathsGenerator(dim=dim, div=div, links=gates, max_cdist=1)
            for paths in paths_gen.generate_paths(std=True):
                fcurve = PathFuzzyCurve.init_from_paths(paths, disable_time_rev=True)
                for curve in fcurve.gen_possible_curves():
                    if all(bm * curve not in seen for bm in all_bms):
                        seen.add(curve)
        self.assertEqual(len(seen), 920)

    def test_junc(self) -> None:
        def is_specialization(curve, tmpl):
            # Check if curve has all of defined specs of the given template, and they are the same."""
            return all(curve.patterns[pnum].specs[cnum] == sp for pnum, cnum, sp in tmpl.gen_defined_specs())

        for num, curve in enumerate(self.curves):
            if num == 0:
                pcurve = curve.forget()
            else:
                pcurve = curve.forget(disable_time_rev=True)

            junc_info = pcurve.get_junction_templates()
            for curve in pcurve.gen_possible_curves():
                juncs = list(curve.gen_regular_junctions())
                # check that for all curve juncs there is a template
                for junc in juncs:
                    self.assertIn(junc, junc_info)
                    self.assertTrue(any(is_specialization(curve, tmpl) for tmpl in junc_info[junc]))

                # check that is curve is consistent with template, then it has junc
                for junc, curves in junc_info.items():
                    if any(is_specialization(curve, tmpl) for tmpl in junc_info[junc]):
                        self.assertIn(junc, juncs)


class TestMisc(unittest.TestCase):
    def test_diag(self) -> None:
        gates = Link.parse_gates('(0,0)->(1,1/2)')
        pgen = PathsGenerator(dim=2, div=5, links=[gates], max_cdist=1)
        paths = next(pgen.generate_paths())
        pcurve = PathFuzzyCurve.init_from_paths(paths)
        curve = pcurve.get_curve_example()
        print(curve)
