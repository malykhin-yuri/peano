"""
This module provides classes for Peano multifractal curves.

Most methods are implemented for the FuzzyCurve class --
curves with some specs undefined. E.g., if specs for
first/last cubes are defined, then you already can
determine entrance/exit points.
"""

from __future__ import annotations
from collections import namedtuple, defaultdict
from dataclasses import dataclass
import itertools
from typing import Iterable, Sequence, Self, Any

from quicktions import Fraction  # type: ignore

from .base_maps import BaseMap, Spec
from .utils import get_periodic_sum
from .paths import Proto, Path, Link, CubeType
from .subsets import Point


class Pattern(namedtuple('_Pattern', ['proto', 'specs'])):
    """
    Pattern - one component of a multifractal, prototype + specifications.
    """

    def __new__(cls, proto: Proto, specs: Iterable[Spec | None]) -> Pattern:
        return super().__new__(cls, proto=proto, specs=tuple(specs))

    @classmethod
    def parse(cls, chain: str, specs: str | Iterable[str]) -> Pattern:
        """
        Parse pattern from bases string.

        Args:
            chain: chain code, see Proto.parse
            specs: list/csv of spec texts, see Spec.parse

        Returns:
            Pattern instance
        """
        proto = Proto.parse(chain)
        specs_strs: list[str]
        if isinstance(specs, str):
            specs_strs = specs.split(',')
        else:
            specs_strs = list(specs)
        specs_list = [Spec.parse(c) for c in specs_strs]
        return cls(proto, specs_list)

    def __str__(self) -> str:
        return '{} | {}'.format(self.proto, ','.join(str(spec) for spec in self.specs))

    def __invert__(self) -> Pattern:
        """Time-reversed pattern."""
        # each base_map does not change:
        #   - pnums are the same
        #   - if there was not time_rev, there should not be after reversal
        #   - cube map does not change, because geometry does not change
        return Pattern(~self.proto, reversed(self.specs))


class FuzzyCurve:
    """
    Multifractal Peano curve that may be not fully specified.

    Peano multifractal consists of several patterns (with the same div!).

    This is abstract class for fuzzy curves.
    The implementation must provide that specs for different fractions
    are independent. This is required for the correct usage of SAT solvers!

    Most methods may raise KeyError if some required specs/proto cubes are not defined
    """

    type PatternArgType = Pattern | tuple[Proto | Iterable[CubeType | None], Iterable[Spec | None]]  # TODO simplify; None in CubeType ?

    def __init__(self, dim: int, div: int, patterns: Iterable[PatternArgType], pnum: int = 0) -> None:
        """
        Create FuzzyCurve instance.

        Args:
            dim: dimension d of cube [0,1]^d (image of the curve)
            div: number of divisions for each of the coordinates, so genus = G = div**dim
            patterns: list of patterns (Pattern instances) or plain (proto, specs) tuples
              some specs in patterns may be None
            pnum: selected pattern, to pick actual curve f:[0,1]->[0,1]^d
        """
        self.dim = dim
        self.div = div

        pattern_objects = []
        for pattern in patterns:
            if not isinstance(pattern, Pattern):
                proto, specs = pattern
                proto = proto if isinstance(proto, Proto) else Proto(dim, div, proto)
                pattern = Pattern(proto=proto, specs=specs)
            pattern_objects.append(pattern)

        self.patterns = tuple(pattern_objects)
        self.mult = len(self.patterns)
        self.pnum = pnum
        self.genus = div**dim

    @property
    def proto(self) -> Proto:
        return self.patterns[self.pnum].proto

    @property
    def specs(self) -> tuple[Spec, ...]:
        return self.patterns[self.pnum].specs

    def __str__(self) -> str:
        return '\n'.join('@{}: {}'.format(pnum, pattern) for pnum, pattern in enumerate(self.patterns))

    def _changed(self, patterns: Iterable[PatternArgType] | None = None, pnum: int | None = None, **kwargs) -> Self:
        return type(self)(
            dim=self.dim,
            div=self.div,
            patterns=patterns if patterns is not None else self.patterns,
            pnum=pnum if pnum is not None else self.pnum,
            **kwargs,
        )

    def __invert__(self) -> Self:
        """Reverse time in a curve."""
        return self._changed(patterns=[~pattern for pattern in self.patterns])

    def _apply_cube_map(self, base_map: BaseMap) -> Self:
        new_patterns: list[FuzzyCurve.PatternArgType] = []
        for pattern in self.patterns:
            # isometry for the prototype
            new_proto = base_map * pattern.proto

            # specs are conjugated: to get from mapped curve its fraction, we can:
            # - map curve to original (base_map^{-1})
            # - then map curve to the fraction (spec)
            # - then map the whole thing to mapped curve (base_map)
            new_specs = [spec.conjugate_by(base_map) if spec is not None else None for spec in pattern.specs]
            new_patterns.append((new_proto, new_specs))

        return self._changed(patterns=new_patterns)

    def __rmul__(self, other: BaseMap | Spec) -> Self:
        """
        Apply base map or spec to a fractal curve, return new curve.
        
        Note that spec also may change curve's selected pnum.
        """
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
        curve = curve._apply_cube_map(base_map.cube_map())

        if curve.pnum != pnum:
            curve = curve._changed(pnum=pnum)

        return curve

    def _compose_specs(self, spec: Spec, cnum: int) -> Spec:
        # Fast spec composition without curve multiplication.
        # Get spec X such that: X * self = C.specs[cnum] * C, where C := spec*self
        # Method allows to get orientations of deep fractions of a curve.
        active_cnum = spec.base_map.apply_cnum(self.genus, cnum)
        last_spec = self.patterns[spec.pnum].specs[active_cnum]
        if last_spec is None:
            raise KeyError
        # pnum is taken from the last spec
        # base_maps: # C.specs[cnum].base_map = spec.base_map * last_spec.base_map * ~spec.base_map, see __rmul__
        #  equivalently: first, we map self by last_spec then apply spec.base_map to the whole picture
        return spec.base_map * last_spec

    def get_deep_spec(self, pnum: int, cnums: Iterable[int]) -> Spec:
        """
        Get spec X such that fraction = X * self.

        Fraction is defined by pnum and nested sequence of cube times (cnums).
        """
        spec = Spec(BaseMap.id_map(self.dim), pnum)
        for cnum in cnums:
            spec = self._compose_specs(spec, cnum)
        return spec

    def gen_allowed_specs(self, pnum, cnum):
        raise NotImplementedError("Define in child class")

    # TODO: return Curve
    def get_curve_example(self) -> Self:
        """
        Generate a curve, compatible with self.

        We use gen_allowed_specs (defined in child classes) to specify missing specs.

        Returns:
            Curve instance
        """
        spec_dict: dict[int, dict[int, Spec]] = defaultdict(dict)
        for pnum, pattern in enumerate(self.patterns):
            for cnum, spec in enumerate(pattern.specs):
                if spec is None:
                    spec_dict[pnum][cnum] = next(self.gen_allowed_specs(pnum=pnum, cnum=cnum))
        return self._specify_allowed(spec_dict)

    def gen_possible_curves(self) -> Iterable[Curve]:
        """Generate all curves, compatible with self."""

        sp_variant_generators = []
        G = self.genus

        # large list of generators of specs for each (pnum, cnum); flatten by pnum * G + cnum
        for pnum in range(self.mult):
            for cnum in range(G):
                sp_variant_generators.append(self.gen_allowed_specs(pnum, cnum))

        for all_specs in itertools.product(*sp_variant_generators):
            patterns = []
            for pnum, pattern in enumerate(self.patterns):
                specs = all_specs[(G * pnum):(G * (pnum + 1))]
                patterns.append((pattern.proto, specs))

            yield Curve(dim=self.dim, div=self.div, patterns=patterns, pnum=self.pnum)

    def count_possible_curves(self) -> int:
        """Count all curves, compatible with self."""
        result = 1
        for pnum in range(self.mult):
            for cnum in range(self.genus):
                result *= len(list(self.gen_allowed_specs(pnum, cnum)))
        return result

    def gen_defined_specs(self) -> Iterable[tuple[int, int, Spec]]:
        """Generate triples (pnum, cnum, spec) of defined specs."""
        for pnum, pattern in enumerate(self.patterns):
            for cnum, spec in enumerate(pattern.specs):
                if spec is not None:
                    yield pnum, cnum, spec

    def specify(self, pnum: int, cnum: int, spec: Spec) -> Self:
        """
        Check that we can set self.patterns[pnum][cnum] = spec, and return specified curve if so.

        Args:
            pnum, cnum: position to specify
            spec: value to specify

        Returns:
            specified curve

        Raises:
            ValueError: if spec is not allowed
        """
        # This is the main method while dividing pairs_tree in estimators,
        # so the efficiency is important here!
        if spec not in self.gen_allowed_specs(pnum, cnum):
            raise ValueError("Can't specify curve")
        if self.patterns[pnum].specs[cnum] is not None:
            return self  # optimization: spec is already defined

        return self._specify_allowed({pnum: {cnum: spec}})

    def _specify_allowed(self, spec_dict: dict[int, dict[int, Spec]]) -> Self:
        # Specify many specs and do not check allowance
        new_patterns = list(self.patterns)
        for pnum, cnum_specs in spec_dict.items():
            pattern = new_patterns[pnum]
            new_specs = list(pattern.specs)
            for cnum, spec in cnum_specs.items():
                new_specs[cnum] = spec
            new_patterns[pnum] = Pattern(pattern.proto, new_specs)
        return self._changed(patterns=new_patterns)

    #
    # Entrance/exit/moments
    #

    def get_entrance(self, pnum: int | None = None) -> Point:
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

    def get_exit(self, pnum: int | None = None) -> Point:
        """Exit of a curve, i.e. point f(1); see get_entrance."""
        if pnum is None:
            pnum = self.pnum
        return self._get_cube_limit(pnum, self.genus-1)

    def _get_cube_limit(self, pnum: int, cnum: int) -> Point:
        # we found the sequence of cubes that we obtain if we take cube #cnum in each fraction
        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=pnum)
        cubes: list[tuple[int, ...]] = []
        index: dict[Spec, int] = {}
        # current curve = cur_spec * self
        while cur_spec not in index:
            index[cur_spec] = len(cubes)
            cube = (cur_spec.base_map * self.patterns[cur_spec.pnum].proto)[cnum]
            if cube is None:
                raise KeyError("Curve not specified enough to get cubes sequence!")
            cubes.append(cube)
            cur_spec = self._compose_specs(cur_spec, cnum)

        idx = index[cur_spec]
        return self._get_limit_point(cubes[0:idx], cubes[idx:])

    def _get_limit_point(self, cubes_start: Sequence[CubeType], cubes_period: Sequence[CubeType]) -> Point:
        xs = []
        for j in range(self.dim):
            xj = get_periodic_sum([cube[j] for cube in cubes_start], [cube[j] for cube in cubes_period], self.div)
            xs.append(xj)
        return Point(xs)

    def get_vertex_moments(self, pnum: int | None = None) -> dict[tuple[int, ...], Fraction]:
        """
        Get all vertex moments, i.e. t: f(t)=v, where v in {0,1}^d

        Note that any vertex moment is uniquely defined
        because only one fraction contains given vertex.

        Args:
            pnum: select non-default pattern

        Returns:
            dict {vertex: moment}.
        """
        return {vertex: self.get_face_moment(vertex, pnum) for vertex in itertools.product((0, 1), repeat=self.dim)}

    def get_face_moment(self, face: tuple[int | None, ...], pnum: int | None = None, last: bool = False, find_point: bool = False) -> Fraction | CurvePoint:
        """
        Moment of face touch (default: first moment).

        Args:
            face: a tuple of {0,1,None} defining the cube face
              {(x_0,...,x_{d-1}): x_i==0 if face[i]==0, x_i==1 if face[i]==1, or arbitrary x[i] if face[i] is None.
              E.g., tuples (0,0,0) or (0,1,1) define vertices
            pnum: select non-default pattern
            last: get last moment instead of first
            find_point: calculate point of face touch

        Returns:
            rational t or CurvePoint instance (if find_point)
        """
        if pnum is None:
            pnum = self.pnum

        cur_spec = Spec(base_map=BaseMap.id_map(self.dim), pnum=pnum)
        cnums: list[int] = []
        cubes: list[CubeType] = []
        index: dict[Spec, int] = {}
        while cur_spec not in index:
            index[cur_spec] = len(cnums)
            cur_proto = cur_spec.base_map * self.patterns[cur_spec.pnum].proto
            cnum = self._get_face_cnum(cur_proto, face, last=last)
            cnums.append(cnum)
            if find_point:
                cubes.append(cur_proto[cnum])
            cur_spec = self._compose_specs(cur_spec, cnum)

        period_start = index[cur_spec]
        moment = get_periodic_sum(cnums[0:period_start], cnums[period_start:], self.genus)
        if find_point:
            point = self._get_limit_point(cubes[0:period_start], cubes[period_start:])
            return CurvePoint(point, moment)
        else:
            return moment

    @staticmethod
    def _get_face_cnum(proto: Proto, face: Sequence[int | None], last: bool = False) -> int:
        # cnum of face touch
        data: Iterable[Any]
        data = enumerate(proto)
        if last:
            data = reversed(list(data))
        for cnum, cube in data:
            touches_face = True
            for x, e in zip(cube, face):
                if (e == 1 and x != (proto.div-1)) or (e == 0 and x != 0):
                    touches_face = False
                    break
            if touches_face:
                return cnum
        raise ValueError("Could not touch face!")  # for mypy :(

    #
    # Junctions; see also the class Junction
    #

    def gen_auto_junctions(self) -> Iterable[AutoJunction]:
        for pnum in range(self.mult):
            yield AutoJunction(dim=self.dim, pnum=pnum)

    def _gen_junctions_from_base(self, base_juncs: Iterable[RegularJunction]) -> Iterable[RegularJunction]:
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

    def _get_base_junction(self, pnum: int, cnum: int, depth: int | None = 1) -> RegularJunction:
        # Get base junction if both specs are defined
        pattern = self.patterns[pnum]
        spec1, spec2 = pattern.specs[cnum], pattern.specs[cnum + 1]
        if spec1 is None or spec2 is None:
            raise KeyError("Can't get base junction: spec is not defined!")
        delta_x = [c2j - c1j for c1j, c2j in zip(pattern.proto[cnum], pattern.proto[cnum + 1])]
        return RegularJunction(spec1, spec2, delta_x, depth=depth)

    def _get_derived_junction(self, junc: RegularJunction) -> RegularJunction:
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

    def get_junction_templates(self) -> dict[RegularJunction, list[FuzzyCurve]]:
        """
        Get possible junctions and template curves for them.

        If a specified curve has some junc, then it is consistent with
        one of the corresponding template curves.
        If a specified curve is consistent with some of the templates,
        then it has the corresponding junc.

        Returns:
            dict {junc: template_curves}
        """
        # template is a tuple (pnum, cnum, tmpl_specs), where in tmpl we specify:
        # - all begin and end specs (to get derivatives)
        # - specs at (pnum, cnum) and (pnum, cnum+1)
        templates = []
        G = self.genus
        P = self.mult
        variants = []
        for pnum in range(P):
            variants.append(self.gen_allowed_specs(pnum=pnum, cnum=0))
            variants.append(self.gen_allowed_specs(pnum=pnum, cnum=G - 1))

        for variant_tuple in itertools.product(*variants):
            variant = list(variant_tuple)
            gate_specs = {}
            for pnum in range(P):
                gate_specs[pnum, 0] = variant.pop(0)
                gate_specs[pnum, G - 1] = variant.pop(0)

            for pnum, cnum in itertools.product(range(P), range(G-1)):
                pc1 = (pnum, cnum)
                pc2 = (pnum, cnum + 1)
                for sp1, sp2 in itertools.product(self.gen_allowed_specs(*pc1), self.gen_allowed_specs(*pc2)):
                    if any(pc in gate_specs and gate_specs[pc] != sp for pc, sp in [(pc1, sp1), (pc2, sp2)]):
                        continue
                    tmpl_specs = gate_specs.copy()
                    tmpl_specs[pc1] = sp1
                    tmpl_specs[pc2] = sp2
                    templates.append((pnum, cnum, tmpl_specs))

        junc_curves = defaultdict(list)
        for pnum, cnum, tmpl_specs in templates:
            spec_dict: dict[int, dict[int, Spec]] = defaultdict(dict)
            for (p, c), sp in tmpl_specs.items():
                spec_dict[p][c] = sp
            curve = self._specify_allowed(spec_dict)
            base_junc = curve._get_base_junction(pnum=pnum, cnum=cnum, depth=None)  # junc depths can't be known for fuzzy curves
            for junc in curve._gen_junctions_from_base([base_junc]):
                junc_curves[junc].append(curve)

        return junc_curves

    def gen_regular_junctions(self) -> Iterable[RegularJunction]:
        """Gen all possible regular junctions."""
        yield from self.get_junction_templates()


class Curve(FuzzyCurve):
    """
    Fully-specified regular Peano curve.

    This curve defines the continuous surjective map f:[0,1]->[0,1]^d.
    """

    def __eq__(self, other):
        return self.patterns, self.pnum == other.patterns, other.pnum

    def __hash__(self):
        return hash((self.patterns, self.pnum))

    @classmethod
    def parse(cls, patterns_bases: Sequence[tuple[str, str | Iterable[str]]]) -> Curve:
        """
        Convenient way to define a curve.

        Args:
            patterns_bases: list of tuples (chain_code, spec_bases) for Pattern.parse

        Returns:
            Curve instance
        """
        patterns = [Pattern.parse(chain, spec_bases) for chain, spec_bases in patterns_bases]
        proto0 = patterns[0].proto
        dim, div = proto0.dim, proto0.div
        return cls(dim, div, patterns)

    def gen_allowed_specs(self, pnum: int, cnum: int) -> Iterable[Spec]:
        yield self.patterns[pnum].specs[cnum]

    def get_paths(self) -> tuple[Path[Point], ...]:
        """Get curve pointed prototypes (Path objects tuple)."""
        links = [Link(self.get_entrance(pnum), self.get_exit(pnum)) for pnum in range(self.mult)]
        paths: list[Path[Point]] = []
        for pattern in self.patterns:
            pattern_links = [spec.base_map * links[spec.pnum] for spec in pattern.specs]
            paths.append(Path(pattern.proto, pattern_links))
        return tuple(paths)

    def forget(self, **kwargs) -> PathFuzzyCurve:
        """Convert curve to a fuzzy curve, saving entrance/exit and forgetting all specs."""
        return PathFuzzyCurve.init_from_paths(self.get_paths(), **kwargs)

    def get_subdivision(self, k: int = 1) -> Curve:
        """Get k-th subdivision of a curve."""
        N = self.div
        current_curve = self
        for _ in range(k):
            new_patterns = []
            for pnum, curr_pattern in enumerate(current_curve.patterns):
                new_proto = []
                new_specs = []
                for cube, spec in zip(curr_pattern.proto, curr_pattern.specs):
                    orig_proto, orig_specs = self.patterns[spec.pnum]  # from original curve
                    new_proto += [tuple(cj*N + ncj for cj, ncj in zip(cube, c)) for c in spec.base_map * orig_proto]
                    if spec.base_map.time_rev:
                        orig_specs = reversed(orig_specs)
                    new_specs += [spec.base_map * sp for sp in orig_specs]
                new_patterns.append((new_proto, new_specs))

            current_curve = Curve(
                dim=self.dim,
                div=current_curve.div * N,  # we change div so do not use ``_changed''
                patterns=new_patterns,
                pnum=self.pnum,
            )

        return current_curve

    def check(self) -> None:
        """Assert that curve (all patterns) is continuous."""
        if any(not path.is_continuous() for path in self.get_paths()):
            raise ValueError("Not contiuous!")

    def _gen_base_junctions(self) -> Iterable[RegularJunction]:
        # junctions from first subdivision
        seen = set()
        for pnum in range(self.mult):
            for cnum in range(self.genus - 1):
                junc = self._get_base_junction(pnum=pnum, cnum=cnum)
                if junc not in seen:
                    yield junc
                    seen.add(junc)

    def gen_regular_junctions(self) -> Iterable[RegularJunction]:
        """Generate all regular junctions for a curve."""
        yield from self._gen_junctions_from_base(self._gen_base_junctions())

    def get_depth(self):  # problems with typing
        return max(junc.depth for junc in self.gen_regular_junctions())


class PathFuzzyCurve(FuzzyCurve):
    """
    Fuzzy curve with fixed paths.

    We fix only paths: prototypes and entrance/exit gates (links) and allow
    the curve to have any specs that are consistent with them.
    Hence the specs choice is independent, as required for applications.
    Some specs may be defined using specify, e.g.

    Additional attributes:
    We store information about "standard" links (gate pairs) used in paths:
    .links_symmetries  --  dict {std_link: symmetries} - list of bms that save std_link
    .links_std  --  dict {std_link: {pnum: std_map}}, such that std_map * pattern_global_link = std_link

    Information about patterns:
    .pattern_links  --  array of std_links for each fraction of each pattern;
                        i.e., if there are two types of links: side (0,0)->(1,0) and diag (0,0)->(1,1),
                        we store here only "type" (side|diag) for each fraction
    .pattern_reprs  --  array of reprs for each pattern; reprs[cnum] = bm that sends std_link to link of fraction
    """

    def __init__(self, *args, **kwargs) -> None:
        self.links_symmetries = kwargs.pop('links_symmetries')
        self.links_std = kwargs.pop('links_std')
        self.pattern_links = kwargs.pop('pattern_links')
        self.pattern_reprs = kwargs.pop('pattern_reprs')
        super().__init__(*args, **kwargs)

    def _changed(self, *args, **kwargs) -> Self:
        for add_field in ['links_symmetries', 'links_std', 'pattern_links', 'pattern_reprs']:
            kwargs.setdefault(add_field, getattr(self, add_field))
        return super()._changed(*args, **kwargs)

    def __invert__(self) -> Self:
        # we do not change links to keep them standard (hence, links_symmetries also do not change)
        new_links_std = {}
        for link, data in self.links_std.items():
            new_links_std[link] = {pnum: ~std_map for pnum, std_map in data.items()}
        new_pattern_links = [list(reversed(links)) for links in self.pattern_links]  # we do not change links, only their order
        new_pattern_reprs = [[~bm for bm in reversed(reprs)] for reprs in self.pattern_reprs]
        # after super invert the curve will be in inconsistent state:(
        return super().__invert__()._changed(
            links_std=new_links_std,
            pattern_links=new_pattern_links,
            pattern_reprs=new_pattern_reprs,
        )

    def _apply_cube_map(self, cube_map: BaseMap) -> Self:
        cube_inv = cube_map**(-1)

        # to get std_link, revert cube_map and the apply old std_map
        new_links_std = {}
        for link, data in self.links_std.items():
            new_links_std[link] = {pnum: std_map * cube_inv for pnum, std_map in data.items()}

        # to get actual link, first map std_link to get old link, the apply cube_map
        new_pattern_reprs = [[cube_map * repr for repr in reprs] for reprs in self.pattern_reprs]

        return super()._apply_cube_map(cube_map)._changed(
            links_std=new_links_std,
            pattern_reprs=new_pattern_reprs,
        )

    def gen_allowed_specs(self, pnum: int, cnum: int) -> Iterable[Spec]:
        pattern = self.patterns[pnum]
        if pattern.specs[cnum] is not None:
            # this is the case for partially specified curves, they arise in dilation estimation
            yield pattern.specs[cnum]
            return

        link = self.pattern_links[pnum][cnum]
        repr = self.pattern_reprs[pnum][cnum]
        symmetries = self.links_symmetries[link]
        std = self.links_std[link]
        for pn, std_map in std.items():
            for symm in symmetries:
                # how do we put pn-pattern in fraction:
                # 1) map pn-pattern link to std_link
                # 2) apply symmetries for std_link
                # 3) apply repr to get fraction link
                yield Spec(repr * symm * std_map, pn)

    @classmethod
    def init_from_paths(cls, paths: Sequence[Path[Point]], base_maps_group: Iterable[BaseMap] | None = None, disable_time_rev: bool = False):
        """
        Create PathFuzzyCurve from a tuple of pointed paths (Path instances).

        Args:
            paths: path for each pattern (so, mult == len(paths))
            base_maps_group: subgroup of all base maps to use in curve
            disable_time_rev: boolean, disable time_rev in curve base_maps
        """

        dim = paths[0].dim
        div = paths[0].div

        if base_maps_group is None:
            base_maps_group = list(BaseMap.gen_base_maps(dim))
        if disable_time_rev:
            base_maps_group = [bm for bm in base_maps_group if not bm.time_rev]

        links_std: dict[Link[Point], dict[int, BaseMap]] = defaultdict(dict)
        for pnum, path in enumerate(paths):
            std_link = min(bm * path.link for bm in base_maps_group)
            std_map = next(bm for bm in base_maps_group if bm * path.link == std_link)
            links_std[std_link][pnum] = std_map

        links_symmetries = {link: [bm for bm in base_maps_group if bm * link == link] for link in links_std}

        patterns, pattern_links, pattern_reprs = [], [], []
        for path in paths:
            specs = (None,) * len(path.proto)  # nothing is defined
            patterns.append((path.proto, specs))

            reprs, links = [], []
            for link in path.links:
                std_link = min(bm * link for bm in base_maps_group)
                if std_link not in links_std:
                    # all links must be allowed images of global (path.link) links!
                    raise KeyError("Found link in path.links that does not correspond to global links!")
                repr0 = next(bm for bm in base_maps_group if bm * std_link == link)
                reprs.append(repr0)
                links.append(std_link)
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
        spec1: first fraction is spec1 * curve
        spec2: second fraction is spec2 * curve
        delta_x: shift vector to get 2-nd fraction from 1-st, element of {0,1,-1}^d
        delta_t: time shift (=0 or 1, see below)
        depth: junction has depth k if it is obtained from fractions of k-th curve subdivision
    """
    def __init__(self, spec1: Spec, spec2: Spec, delta_x: tuple[Fraction, ...], delta_t: int, depth: int | None = None):
        self.spec1 = spec1
        self.spec2 = spec2
        self.delta_x = delta_x
        self.delta_t = delta_t
        self._hash = hash(self._data())  # depth is just some additional information that is determined by junc
        self.depth = depth

    def _data(self) -> tuple[tuple[Fraction, ...], int, Spec, Spec]:
        return self.delta_x, self.delta_t, self.spec1, self.spec2

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
    def __init__(self, spec1: Spec, spec2: Spec, delta_x: Sequence[Fraction], depth: int | None =None):
        if spec1.pnum > spec2.pnum \
                or (spec1.pnum == spec2.pnum and spec1.base_map.time_rev and spec2.base_map.time_rev):
            # swap and reverse time
            delta_x = [-dj for dj in delta_x]
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

    Required for correct dilation estimation
    """
    def __init__(self, dim: int, pnum: int) -> None:
        spec = Spec(base_map=BaseMap.id_map(dim), pnum=pnum)
        super().__init__(spec1=spec, spec2=spec, delta_x=(0,) * dim, delta_t=0, depth=0)


@dataclass(frozen=True, order=True)
class CurvePoint:
    """Pair (x,t) of a Point and time moment; i.e. a point on the graph of a curve, f(t)=x"""
    x: Point
    t: Fraction

    def __rmul__(self, base_map: BaseMap) -> CurvePoint:
        return CurvePoint(base_map * self.x, base_map.apply_t(self.t))

    def scale(self, factor: Fraction) -> CurvePoint:
        return CurvePoint(Point(xj * factor for xj in self.x), (factor**self.x.dim) * self.t)
