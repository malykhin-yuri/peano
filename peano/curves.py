"""
This module provides classes for Peano multifractal curves.

Most methods are implemented for the FuzzyCurve class --
curves with some specs undefined. E.g., if specs for
first/last cubes are defined, then you already can
determine entrance/exit points.
"""

from collections import namedtuple, defaultdict
import itertools

from sympy import Rational

from .base_maps import BaseMap, Spec
from .utils import get_periodic_sum
from .paths import Proto, Path, Link
from .subsets import Point


class Pattern(namedtuple('Pattern', ['proto', 'specs'])):
    """
    Pattern - one component of a multifractal, prototype + specifications.
    """

    @classmethod
    def parse(cls, chain, specs):
        """
        Parse pattern from bases string.

        Args:
            chain: chain code, see Proto.parse
            specs: list/csv of spec texts, see Spec.parse

        Returns:
            Pattern instance
        """
        proto = Proto.parse(chain)
        if isinstance(specs, str):
            specs = specs.split(',')
        specs = [Spec.parse(c) for c in specs]
        return cls(proto, specs)

    def __invert__(self):
        """Time-reversed pattern."""
        # each base_map does not change:
        #   - pnums are the same
        #   - if there was not time_rev, there should not be after reversal
        #   - cube map does not change, because geometry does not change
        return type(self)(~self.proto, tuple(reversed(self.specs)))


class FuzzyCurve:
    """
    Multifractal Peano curve, not fully specified.

    Object attributes:
    .genus -- subj

    Most methods may raise KeyError if some required specs/proto cubes are not defined
    """

    def __init__(self, dim, div, patterns, pnum=0):
        """
        Create FuzzyCurve instance.

        Args:
            dim: dimension d of cube [0,1]^d (image of the curve)
            div: number of divisions for each of the coordinates, so genus = G = div**dim
            patterns -- list of patterns (Pattern instances)
                some specs in patterns may be None
            pnum: selected pattern, to pick actual curve f:[0,1]->[0,1]^d

        Note that we assume that each pattern has the same div.
        """
        self.dim = dim
        self.div = div

        self.patterns = []
        for proto, specs in patterns:
            proto = proto if isinstance(proto, Proto) else Proto(dim, div, proto)
            specs = tuple(sp if isinstance(sp, Spec) else Spec(sp) if sp is not None else None for sp in specs)
            pattern = Pattern(proto=proto, specs=specs)
            self.patterns.append(pattern)

        self.pcount = len(self.patterns)
        self.pnum = pnum
        self.genus = div**dim

    @property
    def proto(self):
        return self.patterns[self.pnum].proto

    @property
    def specs(self):
        return self.patterns[self.pnum].specs

    @classmethod
    def parse(cls, patterns_bases):
        """
        Convenient way to define a curve.

        Args:
            patterns_bases: A list of pattern data, see Pattern.parse

        Returns:
            FuzzyCurve instance
        """
        patterns = [Pattern.parse(chain, spec_bases) for chain, spec_bases in patterns_bases]
        proto0 = patterns[0].proto
        dim, div = proto0.dim, proto0.div
        return cls(dim, div, patterns)

    def __eq__(self, other):
        return self.pcount == other.pcount and all(p1 == p2 for p1, p2 in zip(self.patterns, other.patterns))

    def __str__(self):
        lines = []
        for pnum, pattern in enumerate(self.patterns):
            lines.append('@{}: {} | {}\n'.format(pnum, pattern.proto, ','.join([str(spec) for spec in pattern.specs])))
        return ''.join(lines)

    def changed(self, patterns=None, pnum=None, **kwargs):
        """Create curve with changed parameters."""
        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=patterns if patterns is not None else self.patterns,
            pnum=pnum if pnum is not None else self.pnum,
            **kwargs,
        )

    def __invert__(self):
        """Reverse time in a curve."""
        return self.changed(patterns=[~pattern for pattern in self.patterns])

    def apply_cube_map(self, base_map):
        """Apply cube isometry to a curve."""

        if base_map.dim != self.dim:
            raise Exception("Incompatible base map!")
        if base_map.time_rev:
            raise Exception("Do not use this method with time_rev!")

        new_patterns = []
        for pattern in self.patterns:
            # isometry for the prototype
            new_proto = base_map * pattern.proto

            # specs are conjugated: to get from mapped curve its fraction, we can:
            # - map curve to original (base_map^{-1})
            # - then map curve to the fraction (spec)
            # - then map the whole thing to mapped curve (base_map)
            new_specs = [spec.conjugate_by(base_map) if spec is not None else None for spec in pattern.specs]
            new_patterns.append((new_proto, new_specs))

        return self.changed(patterns=new_patterns)

    def __rmul__(self, other):
        """Apply base map or spec to a fractal curve, return new curve."""
        if isinstance(other, Spec):
            base_map = other.base_map
            pnum = other.pnum
        elif isinstance(other, BaseMap):
            base_map = other
            pnum = self.pnum
        else:
            raise NotImplementedError

        # base_map is the composition of commutating maps: cube_map and time_map
        curve = self
        if base_map.time_rev:
            curve = ~curve
        curve = curve.apply_cube_map(base_map.cube_map())

        if curve.pnum != pnum:
            curve = curve.changed(pnum=pnum)

        return curve

    def compose_specs(self, spec, cnum):
        """
        Fast spec composition without curve multiplication.

        Get spec X = C.specs[cnum] * spec, where C := spec*self,
        i.e. this is spec such that: C.specs[cnum] * C = X * self
        Method allows to get orientations of deep fractions of a curve.

        Args:
            spec: specification defining curve C = spec*self
            cnum: cube index in curve C for next spec

        Returns:
            composed Spec
        """

        active_cnum = spec.base_map.apply_cnum(self.genus, cnum)
        last_spec = self.patterns[spec.pnum].specs[active_cnum]
        if last_spec is None:
            raise KeyError
        # pnum is taken from the last spec
        # base_maps:
        # C.specs[cnum].base_map = spec.base_map * last_spec.base_map * ~spec.base_map, see __rmul__
        return spec.base_map * last_spec

    def gen_allowed_specs(self, pnum, cnum):
        raise NotImplementedError("Define in child class")

    def get_curve_example(self):
        """
        Generate a curve, compatible with self.

        We use gen_allowed_specs (defined in child classes) to specify missing specs.

        Returns:
            Curve instance
        """
        current = self
        for pnum in range(self.pcount):
            for cnum in range(self.genus):
                if current.patterns[pnum].specs[cnum] is None:
                    spec = next(current.gen_allowed_specs(pnum=pnum, cnum=cnum))
                    current = current.specify(pnum=pnum, cnum=cnum, spec=spec)
        return current

    # TODO: function not in prod!
    def gen_possible_curves(self):
        """
        Generate all curves, compatible with self.

        We use gen_allowed_specs (defined in child classes) and suppose
        that specs for different fractions are independent. Required for SAT!!!!!!!
        This is very important condition, which is provided in PathFuzzyCurve class by continuity.
        """

        sp_variant_generators = []
        G = self.genus

        # large list of generators of specs for each (pnum, cnum); flatten by pnum * G + cnum
        for pnum in range(self.pcount):
            for cnum in range(G):
                sp_variant_generators.append(self.gen_allowed_specs(pnum, cnum))

        for all_specs in itertools.product(*sp_variant_generators):
            patterns = []
            for pnum, pattern in enumerate(self.patterns):
                specs = all_specs[(G * pnum):(G * (pnum + 1))]
                patterns.append((pattern.proto, specs))

            yield Curve(dim=self.dim, div=self.div, patterns=patterns, pnum=self.pnum)

    def gen_defined_specs(self):
        """Generate triples (pnum, cnum, spec) of defined specs."""
        for pnum, pattern in enumerate(self.patterns):
            for cnum, spec in enumerate(pattern.specs):
                if spec is not None:
                    yield pnum, cnum, spec

    def specify(self, pnum, cnum, spec):
        """
        Check that we can set spec at pnum, cnum, and return specified curve if so.

        Try to set self.patterns[pnum][cnum] = spec

        Args:
            pnum, cnum: position to specify
            spec: value to specify

        Returns:
            new curve

        Raises:
            ValueError: if spec is not allowed
        """
        # This is the main method while dividing pairs_tree in estimators,
        # so the efficiency is important here!
        if spec not in self.gen_allowed_specs(pnum, cnum):
            raise ValueError("Can't specify curve")

        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            return self  # optimization

        new_specs = list(pattern.specs)
        new_specs[cnum] = spec
        new_pattern = (pattern.proto, new_specs)

        new_patterns = list(self.patterns)
        new_patterns[pnum] = new_pattern

        return self.changed(patterns=new_patterns)

    def specify_dict(self, spec_dict):
        """
        Specify many specs.

        Args:
            spec_dict: dict {(pnum, cnum): spec}

        Returns, Raises: see specify
        """
        curve = self
        for (pnum, cnum), spec in spec_dict.items():
            curve = curve.specify(pnum=pnum, cnum=cnum, spec=spec)
        return curve

    #
    # Entrance/exit/moments
    #

    def get_entrance(self, pnum=None):
        """
        Entrance of a curve, i.e. point f(0).

        Args:
            pnum: if set, find the entrance of pattern pnum (default is to use self.pnum)

        Returns:
            Point instance
        """
        if pnum is None:
            pnum = self.pnum
        return self._get_cube_limit(pnum, 0)

    def get_exit(self, pnum=None):
        """Exit of a curve, i.e. point f(1); see get_entrance."""
        if pnum is None:
            pnum = self.pnum
        return self._get_cube_limit(pnum, self.genus-1)

    def _get_cube_limit(self, pnum, cnum):
        # we found the sequence of cubes that we obtain if we take cube #cnum in each fraction
        # returns pair (non-periodic part, periodic part)
        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=pnum)
        cubes = []
        index = {}

        while True:
            # cur_curve = cur_spec * self
            cur_curve_proto = cur_spec.base_map * self.patterns[cur_spec.pnum].proto
            cube = cur_curve_proto[cnum]
            if cube is None:
                raise KeyError("Curve not specified enough to get cubes sequence!")

            cubes.append(cube)
            index[cur_spec] = len(cubes)-1

            # cur_spec = cur_curve.specs[cnum] * cur_spec
            cur_spec = self.compose_specs(cur_spec, cnum)

            if cur_spec in index:
                idx = index[cur_spec]
                start, period = cubes[0:idx], cubes[idx:]
                break

        pt = []
        for j in range(self.dim):
            start_j = [x[j] for x in start]
            period_j = [x[j] for x in period]
            pt.append(get_periodic_sum(start_j, period_j, self.div))
        return Point(pt)

    def get_vertex_moments(self, pnum=None):
        """
        Get all vertex moments.

        Note that vertex moment, i.e. t: f(t)=v, is uniquely defined
        because only one fraction contains given vertex.

        Args:
            pnum: select non-default pattern

        Returns:
            dict {vertex: moment}.
        """
        return {vertex: self.get_face_moment(vertex, pnum) for vertex in itertools.product((0, 1), repeat=self.dim)}

    def get_face_moment(self, face, pnum=None, last=False):
        """
        Moment of face touch (default: first moment).

        Args:
            face: a tuple of {0,1,None} defining the cube face
              {(x_0,...,x_{d-1}): x_i==0 if face[i]==0, x_i==1 if face[i]==1, or arbitrary x[i] if face[i] is None.
              E.g., tuples (0,0,0) or (0,1,1) define vertices
            pnum: select non-default pattern
            last: get last moment instead of first

        Returns:
            rational number, moment of first touch
        """
        if pnum is None:
            pnum = self.pnum

        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=pnum)
        curve = cur_spec * self  # invariant
        cnums = []
        index = {}
        while True:
            cnum = curve._get_face_cnum(face, last=last)
            cnums.append(cnum)
            index[cur_spec] = len(cnums)-1

            sp = curve.specs[cnum]
            cur_spec = sp * cur_spec
            curve = sp * curve

            if cur_spec in index:
                period_start = index[cur_spec]
                return get_periodic_sum(cnums[0:period_start], cnums[period_start:], self.genus)

    def _get_face_cnum(self, face, last=False):
        data = enumerate(self.proto)
        if last:
            data = reversed(list(data))
        for cnum, cube in data:
            # check that cube touches face
            touch = True
            for x, e in zip(cube, face):
                if (e == 1 and x != (self.div-1)) or (e == 0 and x != 0):
                    touch = False
                    break
            if touch:
                return cnum

    #
    # Junctions; see also the class Junction
    #

    def gen_auto_junctions(self):
        return (AutoJunction(dim=self.dim, pnum=pnum) for pnum in range(self.pcount))

    def _gen_junctions_from_base(self, base_juncs):
        # Yield base junctions and their derivatives
        # provide correct derivation order (in-width) to get correct depth in Curve.gen_regular_junctions
        seen = set()
        to_derive = []
        for junc in base_juncs:
            yield junc
            seen.add(junc)
            to_derive.append(junc)

        while to_derive:
            # to_derive is a queue (fifo), this guarantees that depths are correct
            junc = to_derive.pop()
            dj = self._get_derived_junction(junc)
            if dj not in seen:
                yield dj
                seen.add(dj)
                to_derive.append(dj)

    def _get_base_junction(self, pnum, cnum):
        # Get base junction if both specs are defined
        pattern = self.patterns[pnum]
        spec1, spec2 = pattern.specs[cnum], pattern.specs[cnum + 1]
        if spec1 is None or spec2 is None:
            raise KeyError("Can't get base junction: spec is not defined!")
        delta_x = [c2j - c1j for c1j, c2j in zip(pattern.proto[cnum], pattern.proto[cnum + 1])]
        return RegularJunction(spec1, spec2, delta_x, depth=1)

    def _get_derived_junction(self, junc):
        if not isinstance(junc, RegularJunction):
            raise ValueError("Derivative is defined for regular junctions!")

        spec1 = junc.spec1
        spec2 = junc.spec2

        p1 = self.patterns[spec1.pnum]
        p2 = self.patterns[spec2.pnum]

        cnum1 = 0 if spec1.base_map.time_rev else -1  # cnum in self
        cnum2 = -1 if spec2.base_map.time_rev else 0

        if p1.specs[cnum1] is None or p2.specs[cnum2] is None:
            raise KeyError("Can't get derivative: spec not defined")

        cube1 = spec1.base_map.apply_cube(self.div, p1.proto[cnum1])  # spec1.base_map == id, so no need in this
        cube2 = spec2.base_map.apply_cube(self.div, p2.proto[cnum2])
        der_delta = tuple(c2j + dxj * self.div - c1j for c1j, c2j, dxj in zip(cube1, cube2, junc.delta_x))

        return RegularJunction(
            spec1.base_map * p1.specs[cnum1],
            spec2.base_map * p2.specs[cnum2],
            der_delta,
            depth=(junc.depth + 1 if junc.depth is not None else None),
        )

    def get_junctions_info(self):
        """
        Get possible junctions and specs for them.

        Returns:
            dict {junc: curves}, with partially specified curves that have given junction.
        """

        # config is a tuple (pnum, cnum, curve), where in curve there are specified:
        # - all begin and end specs (to get derivatives)
        # - specs at (pnum, cnum) and (pnum, cnum+1)
        configs = []
        G = self.genus
        P = self.pcount
        variants = []
        for pnum in range(P):
            variants.append(self.gen_allowed_specs(pnum=pnum, cnum=0))
            variants.append(self.gen_allowed_specs(pnum=pnum, cnum=G - 1))

        for variant in itertools.product(*variants):
            variant = list(variant)
            gate_specs = {}
            for pnum in range(P):
                gate_specs[(pnum, 0)] = variant.pop(0)
                gate_specs[(pnum, G - 1)] = variant.pop(0)

            for pnum in range(P):
                for cnum in range(G - 1):
                    pc1 = (pnum, cnum)
                    for sp1 in self.gen_allowed_specs(*pc1):
                        if pc1 in gate_specs and gate_specs[pc1] != sp1:
                            continue
                        pc2 = (pnum, cnum + 1)
                        for sp2 in self.gen_allowed_specs(*pc2):
                            if pc2 in gate_specs and gate_specs[pc2] != sp2:
                                continue
                            def_specs = gate_specs.copy()
                            def_specs[pc1] = sp1
                            def_specs[pc2] = sp2
                            configs.append((pnum, cnum, self.specify_dict(def_specs)))

        junc_curves = defaultdict(list)
        for pnum, cnum, curve in configs:
            base_junc = curve._get_base_junction(pnum=pnum, cnum=cnum)
            for junc in curve._gen_junctions_from_base([base_junc]):
                # junc depths may be incorrect for fuzzy curves, because depth depends on all of curve's juncs
                junc.depth = None
                junc_curves[junc].append(curve)

        return junc_curves

    def gen_regular_junctions(self):
        """Gen all possible(!) regular junctions."""
        yield from self.get_junctions_info().keys()


class Curve(FuzzyCurve):
    """
    Fully-specified regular Peano curve.

    This curve defines the continuous surjective map f:[0,1]->[0,1]^d.
    """

    def _gen_base_junctions(self):
        # junctions from first subdivision
        seen = set()
        for pnum in range(self.pcount):
            for cnum in range(self.genus - 1):
                junc = self._get_base_junction(cnum=cnum, pnum=pnum)
                if junc not in seen:
                    yield junc
                    seen.add(junc)

    def gen_regular_junctions(self):
        """Generate all regular junctions for a curve."""
        yield from self._gen_junctions_from_base(self._gen_base_junctions())

    def get_depth(self):
        return max(junc.depth for junc in self.gen_regular_junctions())

    def gen_allowed_specs(self, pnum, cnum):
        yield self.patterns[pnum].specs[cnum]

    def get_paths(self):
        # TODO: unite with check
        gate_pairs = [Link(self.get_entrance(pnum), self.get_exit(pnum)) for pnum in range(self.pcount)]
        paths = []
        for pnum, pattern in enumerate(self.patterns):
            pattern_gate_pairs = [spec.base_map * gate_pairs[spec.pnum] for spec in pattern.specs]
            paths.append(Path(pattern.proto, pattern_gate_pairs))
        return paths

    def forget(self, allow_time_rev=False):
        """Convert curve to a fuzzy curve, saving entrance/exit and forgetting all specs."""
        return PathFuzzyCurve.init_from_paths(self.get_paths(), allow_time_rev=allow_time_rev)

    def get_subdivision(self, k=1):
        """Get k-th subdivision of a curve."""
        N = self.div
        current_curve = self
        for _ in range(k):
            new_patterns = []
            for pnum, curr_pattern in enumerate(current_curve.patterns):
                new_proto = []
                new_specs = []
                for cube, spec in zip(curr_pattern.proto, curr_pattern.specs):
                    proto, specs = self.patterns[spec.pnum]  # from original curve

                    if spec.base_map.time_rev:
                        proto = reversed(proto)
                        specs = reversed(specs)

                    for c in proto:
                        nc = spec.base_map.apply_cube(N, c)
                        new_cube = [cj*N + ncj for cj, ncj in zip(cube, nc)]
                        new_proto.append(new_cube)

                    # базовые преобразования для подраздедения:
                    # пусть (cube, spec) соответствуют i-й фракции
                    # в ней мы взяли j-ю подфракцию (sp)
                    # Какое преобразование переводит кривую в j-ю фракцию внутри i-й?
                    # - сначала к исходной кривой мы применим bm, чтобы перевести её в j-ю фракцию,
                    # - потом ко всей этой картинке применяем base_map, чтобы перевести всё в i-ю фракцию (base_map)
                    # можно сделать наоборот:
                    # - сначала кривую переводим в i-ю фракцию (base_map)
                    # - применяем внутри i-й фракции преобразования для перехода в j-ю
                    #   но там оно будет сопряженное: base_map * bm * base_map^{-1}, см. apply_base_map
                    for sp in specs:
                        new_specs.append(spec.base_map * sp)

                new_patterns.append((new_proto, new_specs))

            current_curve = type(self)(
                dim=self.dim,
                div=N*current_curve.div,  # we change div so do not use ``changed''
                patterns=new_patterns,
            )

        return current_curve

    def check(self):
        # TODO а нужен ли? или в тесты его??
        """
        Check consistency of curve params.

        Main check is the continuity of the curve.
        It is equivalent to exit/entrance correspondence for all of the patterns
        (there is a subtlety for non-active patterns - we check them, too).
        """

        d, n, G, P = self.dim, self.div, self.genus, self.pcount

        for pattern in self.patterns:
            assert len(pattern.proto) == self.genus, 'bad proto length'
            for cube in pattern.proto:
                for j in range(d):
                    assert 0 <= cube[j] < n, 'bad cube coordinates'
            assert len(set(pattern.proto)) == len(pattern.proto), 'non-unique cubes'

        # main check - continuity - for all patterns
        entr = {pnum: self.get_entrance(pnum) for pnum in range(P)}
        exit = {pnum: self.get_exit(pnum) for pnum in range(P)}
        for pnum, pattern in enumerate(self.patterns):
            curve_entr = entr[pnum]
            curve_exit = exit[pnum]

            gates = [(None, curve_entr)]  # pairs (entr, exit) in [0,1]^d

            for cube, spec in zip(pattern.proto, pattern.specs):
                bm = spec.base_map
                frac_gates_rel = [bm * point for point in [entr[spec.pnum], exit[spec.pnum]]]
                frac_gates = [tuple((Rational(c, 1) + e) * Rational(1, n) for c, e in zip(cube, point)) for point in frac_gates_rel]
                if bm.time_rev:
                    frac_gates.reverse()
                gates.append(frac_gates)

            gates.append((curve_exit, None))

            for i in range(len(gates)-1):
                if gates[i][1] != gates[i+1][0]:
                    msg = 'exit does not correspond to entrance at ' + str(i)
                    print(gates)
                    raise Exception(msg)


class PathFuzzyCurve(FuzzyCurve):
    """
    Fuzzy curve with fixed paths.

    Additional attributes:
    Information about "standard" links used in paths:
    .links_symmetries  --  dict {std_link: symmetries} - list of bms that save std_link (doest not change at *)
    .links_std  --  dict {std_link: {pnum: std_map}}, such that std_map * pattern_global_link = std_link

    Information about patterns:
    .pattern_links  --  array of std_links for each fraction of each pattern;
    .pattern_reprs  --  array of reprs for each pattern; reprs[cnum] = bm that sends std_link to link of fraction
    """

    def __init__(self, *args, **kwargs):
        links_symmetries = kwargs.pop('links_symmetries')
        links_std = kwargs.pop('links_std')
        pattern_links = kwargs.pop('pattern_links')
        pattern_reprs = kwargs.pop('pattern_reprs')
        super().__init__(*args, **kwargs)
        self.links_symmetries = links_symmetries
        self.links_std = links_std
        self.pattern_links = pattern_links
        self.pattern_reprs = pattern_reprs

    def changed(self, *args, **kwargs):
        for add_field in ['links_symmetries', 'links_std', 'pattern_links', 'pattern_reprs']:
            if add_field not in kwargs:
                kwargs[add_field] = getattr(self, add_field)
        return super().changed(*args, **kwargs)

    def __invert__(self):
        # we do not change links to keep them standard (symmetries also do not change)
        new_links_std = {}
        for link, data in self.links_std.items():
            new_links_std[link] = {pnum: ~std_map for pnum, std_map in data.items()}

        new_pattern_links = [list(reversed(links)) for links in self.pattern_links]
        new_pattern_reprs = []
        for reprs in self.pattern_reprs:
            new_reprs = list(reversed([~bm for bm in reprs]))
            new_pattern_reprs.append(new_reprs)
        return super().__invert__().changed(
            links_std=new_links_std,
            pattern_links=new_pattern_links,
            pattern_reprs=new_pattern_reprs,
        )

    def apply_cube_map(self, cube_map):
        curve = super().apply_cube_map(cube_map)
        new_links_std = {}
        for link, data in self.links_std.items():
            new_links_std[link] = {pnum: std_map * cube_map**(-1) for pnum, std_map in data.items()}

        new_pattern_reprs = []
        for reprs in self.pattern_reprs:
            new_reprs = [cube_map * repr for repr in reprs]
            new_pattern_reprs.append(new_reprs)

        return curve.changed(
            links_std=new_links_std,
            pattern_reprs=new_pattern_reprs,
        )

    def gen_allowed_specs(self, pnum, cnum):
        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            yield pattern.specs[cnum]
            return

        link = self.pattern_links[pnum][cnum]
        repr = self.pattern_reprs[pnum][cnum]
        symmetries = self.links_symmetries[link]
        std = self.links_std[link]
        for pnum in sorted(std.keys()):
            for symm in symmetries:
                # how go we get fraction from a pattern:
                # 1) map pattern link to std_link
                # 2) apply symmetries for std_link
                # 3) apply repr to get fraction link
                yield Spec(repr * symm * std[pnum], pnum)

    @classmethod
    def init_from_paths(cls, paths, allow_time_rev=True):
        dim = paths[0].dim
        div = paths[0].div
        if allow_time_rev:
            possible_maps = list(BaseMap.gen_base_maps(dim))
        else:
            possible_maps = list(BaseMap.gen_base_maps(dim, time_rev=False))
        links_symmetries = {}
        links_std = {}
        for pnum, path in enumerate(paths):
            # pnum stands for path_num and also pattern_num
            std_link = path.link.std()
            std_map = next(bm for bm in possible_maps if bm * path.link == std_link)
            links_std.setdefault(std_link, {})[pnum] = std_map

        for link in links_std:
            links_symmetries[link] = [bm for bm in possible_maps if bm * link == link]

        patterns = []
        pattern_links = []
        pattern_reprs = []
        for pnum, path in enumerate(paths):
            proto = path.proto
            specs = [None] * len(proto)  # nothing is defined
            patterns.append((proto, specs))

            reprs = []
            links = []
            for link in path.links:
                repr0, link0 = None, None
                for std_link in links_std:
                    allowed = [bm for bm in possible_maps if bm * std_link == link]
                    if allowed:
                        # should be exactly once!
                        repr0 = allowed[0]
                        link0 = std_link
                        break
                if repr0 is None or link0 is None:
                    raise Exception("Can't create PathFuzzyCurve: no allowed links")
                reprs.append(repr0)
                links.append(link0)
            pattern_reprs.append(reprs)
            pattern_links.append(links)

        return cls(
            dim=dim, div=div, patterns=patterns,
            links_symmetries=links_symmetries,
            links_std=links_std,
            pattern_links=pattern_links,
            pattern_reprs=pattern_reprs,
        )


class Junction:
    """
    Junction is a pair of two curve fractions.

    Attributes:
    .spec1:  first fraction is spec1 * curve
    .spec2:  second fractions is spec2 * curve
    .delta_x:  shift vector to get 2-nd fraction from 1-st, element of {0,1,-1}^d
    .delta_t:  time shift (=0 or 1, see below)
    .depth:  junction has depth k if it is obtained from fractions of k-th curve subdivision
    """
    def __init__(self, spec1, spec2, delta_x, delta_t, depth=None):
        self.spec1 = spec1
        self.spec2 = spec2
        self.delta_x = delta_x
        self.delta_t = delta_t
        self._hash = hash(self._data())  # depth is just some additional information that is determined by junc
        self.depth = depth

    def _data(self):
        return self.delta_x, self.spec1, self.spec2

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return self._hash

    def __str__(self):
        return '{} | dx={}, dt={} | {}'.format(self.spec1, self.delta_x, self.delta_t, self.spec2)


class RegularJunction(Junction):
    """
    Pair of adjacent curve fractions.

    Junction is standartized:
    - pnum1 <= pnum2
    - first spec has cube_map = id (but possibly with time_rev - this allows to get delta_t=1)
    - always delta_t = 1

    Sources of regular junctions:
    * base junctions
    * derivatives of base junctions (derivative is implemented in the curve class)
    """
    def __init__(self, spec1, spec2, delta_x, depth=None):
        if spec1.pnum > spec2.pnum \
                or (spec1.pnum == spec2.pnum and spec1.base_map.time_rev and spec2.base_map.time_rev):
            # swap and reverse time
            delta_x = tuple(-dj for dj in delta_x)
            spec1, spec2 = ~spec2, ~spec1

        bm1_cube_inv = spec1.base_map.cube_map()**(-1)
        super().__init__(
            spec1=bm1_cube_inv * spec1,  # only possible time_rev
            spec2=bm1_cube_inv * spec2,
            delta_x=bm1_cube_inv.apply_vec(delta_x),
            delta_t=1,
            depth=depth,
        )


class AutoJunction(Junction):
    """
    Auto junction: fraction with itself.

    delta_x = 0, delta_t = 0

    Required for correct ratio estimation
    """
    def __init__(self, dim, pnum=0):
        spec = Spec(base_map=BaseMap.id_map(dim), pnum=pnum)
        super().__init__(spec1=spec, spec2=spec, delta_x=(0,) * dim, delta_t=0, depth=0)
