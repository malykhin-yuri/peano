import itertools
from collections import namedtuple

from quicktions import Fraction

from .base_maps import BaseMap


class Subset:
    """
    Abstract class for subsets in [0,1]^d (or in R^d)

    Defines interface for usage of subsets in Link class.
    """

    def intersects(self, other):
        raise NotImplementedError

    def map_to_cube(self, div, cube):
        """
        Put self in given cube.

        Return f_c(S), where f_c maps [0,1]^d to cube/div.
        """
        raise NotImplementedError

    def gen_neighbours(self):
        """
        Generate neighbour cubes with set portions.

        Given S, yields pairs of (cube_start, S'), where (nonempty) S' is the portion of S in that cube,
        i.e., (S \cap cube) = S' + cube_start;  cube_start != 0
        """
        raise NotImplementedError

    def divide(self, div):
        """
        Get fractions of this set.

        Given S, return pairs of (cube, S'), where (nonempty) S' is the fraction of S in the cube,
        i.e., (S \cap cube/div) = (S'+cube_start)/div
        """
        raise NotImplementedError

    def __rmul__(self, base_map):
        """Apply cube map to self."""
        raise NotImplementedError

    def std(self):
        """Standartization."""
        return min(bm * self for bm in BaseMap.gen_base_maps(self.dim, time_rev=False))

    def argmul_intersect(self, other, bms=None):
        """Yield bms that bm * self intersects other."""
        if bms is None:
            bms = BaseMap.gen_base_maps(self.dim, time_rev=False)
        for bm in bms:
            if (bm * self).intersects(other):
                yield bm


class Point(tuple, Subset):
    """
    Rational point, implemented as tuple of Fractions.

    Used in Gates (point links) and for curve entrance/exit.
    """
    @property
    def dim(self):
        return len(self)

    def transform(self, scale=None, shift=None):
        """M -> M*scale + shift"""
        if scale is not None:
            scale = Fraction(scale)
        new_pt = []
        for j, pj in enumerate(self):
            if scale is not None:
                pj *= scale
            if shift is not None:
                pj += Fraction(shift[j])
            new_pt.append(pj)
        return Point(new_pt)

    def map_to_cube(self, div, cube):
        return self.transform(shift=cube).transform(scale=Fraction(1, div))

    def intersects(self, other):
        return self == other

    def gen_neighbours(self):
        for cube in self._gen_integer_cubes():
            if all(cj == 0 for cj in cube):
                continue
            yield cube, self.transform(shift=[-cj for cj in cube])

    def divide(self, div):
        scaled = self.transform(scale=div)
        for cube in scaled._gen_integer_cubes():
            if any(cj < 0 or cj >= div for cj in cube):
                continue
            yield cube, scaled.transform(shift=[-cj for cj in cube])

    def _gen_integer_cubes(self):
        # cubes that contain point
        rngs = []
        for pj in self:
            r = pj.numerator // pj.denominator
            rng = (r-1, r) if pj.numerator % pj.denominator == 0 else (r,)
            rngs.append(rng)
        yield from itertools.product(*rngs)

    def __rmul__(self, bm):
        new_x = bm.apply_x(self)
        return Point(new_x)

    def __str__(self):
        return '(' + ','.join(str(pj) for pj in self) + ')'

    @classmethod
    def parse(cls, text):
        """Parse from usual representation, e.g. (0,1/2,-1/3)."""
        pt = [Fraction(token) for token in text.strip().strip('()').split(',')]
        return cls(pt)

    def std(self):
        half = Fraction(1, 2)
        one = Fraction(1)
        new = [pj if pj <= half else one - pj for pj in self]
        new.sort()
        return Point(new)

    def argstd(self):
        yield from BaseMap.gen_constraint_fast(self, self.std())

    def face_dim(self):
        return sum(pj != Fraction(0) and pj != Fraction(1) for pj in self)


class FacetDivSubset(Subset):
    """
    Nested-cubes subset of [0,1]^d cube facet (hyperface)

    This class is used in Links for facet-gated curves search.
    An important property of this family: any two sets are subsets or disjoint.
    The (int)boundary is not included in the set (see below).

    Examples:
        FacetDivSubset(dim=3, div=2, facet=(0, 1))  -  x0=1, 0<x1<1, 0<x2<1
        FacetDivSubset(dim=3, div=2, facet=(0, 1), cubes=[(1,1)])  -  x0=1, 1/2<x1<1, 1/2<x2<1
        FacetDivSubset(dim=3, div=3, facet=(2, 1), cubes=[(0,1)])  -  x2=1, 0<x0<1/3, 1/3<x1<2/3
        FacetDivSubset(dim=2, div=2, facet=(0, 0), cubes=[(0,),(1,)])  -  x0=0, 1/4<x1<1/2
    """

    def __init__(self, dim, div, facet, cubes=()):
        """
        Args:
            dim, div: as usual
            facet: pair (coord, value) defining the cube facet: x_{coord}=value
            cubes: sequence of div-cubes defining a position in the facet
        """
        self.dim = dim
        self.div = div
        self.facet = tuple(facet)
        self.cubes = tuple(tuple(cube) for cube in cubes)

    def is_subset(self, other):
        return self.intersects(other) and len(self.cubes) >= len(other.cubes)

    def intersects(self, other):
        if self.facet != other.facet:
            return False
        if any(c1 != c2 for c1, c2 in zip(self.cubes, other.cubes)):  # zip shortest
            return False
        return True

    @property
    def facet_coord(self):
        return self.facet[0]

    @property
    def facet_value(self):
        return self.facet[1]

    @property
    def _data(self):
        return self.dim, self.div, self.facet, self.cubes

    def __eq__(self, other):
        return self._data == other._data

    def __hash__(self):
        return hash(self._data)

    def __lt__(self, other):
        return self._data < other._data

    def _map_facet(self, cube_map):
        # bm acts on cube facet
        for j, (k, b) in enumerate(cube_map.coords):  # y[j] = x[k]^b
            if k == self.facet_coord:
                return j, 1 - self.facet_value if b else self.facet_value

    def _abs_to_facet(self, j):
        return j if j < self.facet_coord else j-1

    def __rmul__(self, cube_map):
        """Apply cube map (a BaseMap instance)."""
        assert self.dim == cube_map.dim

        new_facet = self._map_facet(cube_map)

        # j = absolute coordinate, m = relative in self.facet, i = relative in new facet
        # in-facet coordinate map: (u_0,...,u_{facet_dim-1}) -> (v_0,...,v_{facet_dim-1})
        facet_map_data = []
        for j in range(self.dim):
            if j == new_facet[0]:
                continue
            k, b = cube_map.coords[j]  # y_j = x_k^b
            m = self._abs_to_facet(k)  # x_k = u_m
            facet_map_data.append((m, b))  # v_i = y_j = u_m^b  (i=index)

        facet_map = BaseMap(facet_map_data)
        new_cubes = (facet_map.apply_cube(self.div, cube) for cube in self.cubes)
        return FacetDivSubset(self.dim, self.div, new_facet, new_cubes)

    def all_cube_maps(self):
        if not hasattr(self, '_all_cube_maps'):
            self._all_cube_maps = list(BaseMap.gen_base_maps(self.dim, time_rev=False))
        return self._all_cube_maps

    def argmul_intersect(self, other, bms=None):
        """Yield bms that bm * self intersects other."""
        if bms is None:
            bms = self.all_cube_maps()
        for bm in bms:
            if self._map_facet(bm) != other.facet:  # optimization!
                continue
            if (bm * self).intersects(other):
                yield bm

    def map_to_cube(self, div, cube):
        assert div == self.div
        assert cube[self.facet_coord] == (0 if self.facet_value == 0 else div-1)
        facet_cube = tuple(cube[c] for c in range(self.dim) if c != self.facet_coord)
        return FacetDivSubset(self.dim, self.div, self.facet, (facet_cube,) + self.cubes)

    def gen_neighbours(self):
        if self.facet_value == 0:
            delta = -1
            new_facet_value = 1
        else:
            delta = 1
            new_facet_value = 0
        cube = tuple(0 if j != self.facet_coord else delta for j in range(self.dim))
        yield cube, FacetDivSubset(self.dim, self.div, (self.facet_coord, new_facet_value), self.cubes)

    def divide(self, div):
        assert div == self.div
        dim = self.dim
        fc, fv = self.facet

        # facet_cube = (dim-1)-cube on the facet with relative coordinates
        if self.cubes:
            facet_cubes = [self.cubes[0]]
            new_set = FacetDivSubset(dim, div, self.facet, self.cubes[1:])
        else:
            facet_cubes = itertools.product(range(div), repeat=dim-1)
            new_set = self

        for facet_cube in facet_cubes:
            abs_cube = tuple((0 if fv == 0 else div-1) if j == fc else facet_cube[self._abs_to_facet(j)] for j in range(dim))
            yield abs_cube, new_set

    def __str__(self):
        text = 'x{}={}'.format(self.facet_coord, self.facet_value)
        if self.cubes:
            cubes_str = '->'.join(str(cube) for cube in self.cubes)
            text += ':{}'.format(cubes_str)
        return text


class Link(namedtuple('Link', ['entrance', 'exit'])):
    """
    Link is defined by pair of subsets of [0,1]^d: entrance and exit subsets (derived from Subset class)

    Link defines the set of curves with curve.entrance in link.entrance and curve.exit in link.exit.
    Instances of this class may be used in path generators.

    Most useful are pointed links - pairs of entrance/exit gates.
    """
    @property
    def dim(self):
        return self.entrance.dim

    def __rmul__(self, base_map):
        ne, nx = base_map * self.entrance, base_map * self.exit
        if base_map.time_rev:
            ne, nx = nx, ne
        return type(self)(ne, nx)

    def transform(self, *args, **kwargs):
        return type(self)(
            self.entrance.transform(*args, **kwargs),
            self.exit.transform(*args, **kwargs),
        )

    def std(self):
        # heavily-used function, optimize it for point gates
        if not hasattr(self.entrance, 'argstd'):
            return min(bm * self for bm in BaseMap.gen_base_maps(self.dim))

        entr, exit = self.entrance, self.exit
        std1 = entr.std()
        std2 = exit.std()
        if std1 < std2:
            min_std = std1
            other = min(bm * exit for bm in entr.argstd())
        elif std1 > std2:
            min_std = std2
            other = min(bm * entr for bm in exit.argstd())
        else:
            min_std = std1
            other1 = min(bm * exit for bm in entr.argstd())
            other2 = min(bm * entr for bm in exit.argstd())
            other = min(other1, other2)
        return type(self)(min_std, other)

    def __invert__(self):
        return type(self)(self.exit, self.entrance)

    def __str__(self):
        return '{} -> {}'.format(self.entrance, self.exit)

    def intersects(self, other):
        return self.entrance.intersects(other.entrance) and self.exit.intersects(other.exit)

    def argmul_intersect(self, other):
        bms1 = list(self.entrance.argmul_intersect(other.entrance))
        for bm in self.exit.argmul_intersect(other.exit, bms=bms1):
            yield bm
        bms2 = list(self.entrance.argmul_intersect(other.exit))
        for bm in self.exit.argmul_intersect(other.entrance, bms=bms2):
            yield ~bm

    def is_pointed(self):
        """Pair of Points, e.g., entrance and exit gates of a curve."""
        return isinstance(self.entrance, Point) and isinstance(self.exit, Point)

    @classmethod
    def parse_gates(cls, text):
        """Create Points Link from string, e.g., '(0,0)->(1,1/2)'"""
        points = [Point.parse(pt_str) for pt_str in text.split('->')]
        assert len(points) == 2
        return cls(*points)
