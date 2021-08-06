import itertools
from collections import namedtuple

from sympy import Rational

from .base_maps import BaseMap


class Subset:
    """
    Abstract class for subsets in [0,1]^d (or in R^d)

    Defines interface for usage of subsets in Link class.
    """

    def intersects(self, other):
        raise NotImplementedError

    def map_to_cube(self, div, cube):
        """Put self in given cube.

        Return f_c(S), where f_c maps [0,1]^d to cube/div.
        """
        raise NotImplementedError

    def gen_neighbours(self):
        """Generate neighbour cubes with set portions.

        Given S, yields pairs of (cube_start, S'), where (nonempty) S' is the portion of S in that cube,
        i.e., (S \cap cube) = S' + cube_start
        """
        raise NotImplementedError

    def divide(self, div):
        """Get fractions of this set.

        Given S, return pairs of (cube, S'), where (nonempty) S' is the fraction of S in the cube,
        i.e., (S \cap cube/div) = (S'+cube_start)/div
        """
        raise NotImplementedError

    def __rmul__(self, base_map):
        """Apply cube map to self."""
        raise NotImplementedError

    def std(self):
        """Standartization."""
        return min(bm * self for bm in BaseMap.gen_base_maps(self.dim))

    def argmul_intersect(self, other, bms=None):
        """Yield bms that bm * self intersects other."""
        if bms is None:
            bms = BaseMap.gen_base_maps(self.dim, time_rev=False)
        for bm in bms:
            if (bm * self).intersects(other):
                yield bm


class Point(tuple, Subset):
    """Tuple of FastFractions."""

    @property
    def dim(self):
        return len(self)

    def transform(self, scale=None, shift=None):
        """M -> M*scale + shift"""
        if scale is not None:
            scale = Rational(scale)
        new_pt = []
        for j, pj in enumerate(self):
            if scale is not None:
                pj *= scale
            if shift is not None:
                pj += Rational(shift[j])
            new_pt.append(pj)
        return type(self)(new_pt)

    def map_to_cube(self, div, cube):
        return self.transform(shift=cube).transform(scale=Rational(1, div))

    def intersects(self, other):
        return self == other

    def gen_neighbours(self):
        for cube in self.gen_integer_cubes():
            if all(cj == 0 for cj in cube):
                continue
            yield cube, self.transform(shift=[-cj for cj in cube])

    def divide(self, div):
        scaled = self.transform(scale=div)
        for cube in scaled.gen_integer_cubes():
            if any(cj < 0 for cj in cube):
                continue
            if any(cj >= div for cj in cube):
                continue
            yield cube, scaled.transform(shift=[-cj for cj in cube])

    def gen_integer_cubes(self):
        rngs = []
        for pj in self:
            r = pj.p // pj.q
            rng = [r-1, r] if pj.p % pj.q == 0 else [r]
            rngs.append(rng)
        yield from itertools.product(*rngs)

    def __rmul__(self, bm):
        return type(self)(Rational(1) - self[k] if b else self[k] for k, b in bm.coords)

    def __str__(self):
        return '(' + ','.join([str(pj) for pj in self]) + ')'

    @classmethod
    def parse(cls, text):
        pt = [Rational(token) for token in text.strip().strip('()').split(',')]
        return cls(pt)

    def std(self):
        Half = Rational(1, 2)
        new = [pj if pj <= Half else Rational(1) - pj for pj in self]
        new.sort()
        return type(self)(new)

    def argstd(self):
        return list(BaseMap.gen_constraint_fast(self, self.std()))

    def face_dim(self):
        return sum(pj != Rational(0) and pj != Rational(1) for pj in self)


class FacetDivSubset(Subset):
    """
    Nested-cubes subset of [0,1]^d cube facet (hyperface)

    Instances of this class are intended to be used in Portals for hypercurves search.
    An important property of this family: any two sets are subsets or disjoint.

    Examples:
        FacetDivSubset(dim=3, div=2, face=(0: 1))  -  x0=1, 0<x1<1, 0<x2<1
        FacetDivSubset(dim=3, div=2, face=(0: 1), cubes=[(1,1)])  -  x0=1, 1/2<x1<1, 1/2<x2<1
        FacetDivSubset(dim=3, div=3, face=(1: 1), cubes=[(0,1)])  -  x1=1, 0<x0<1/3, 1/3<x2<2/3
        FacetDivSubset(dim=2, div=2, face={0: 0), cubes=[(0,),(1,)])  -  x0=0, 1/4<x1<1/2
    """

    def __init__(self, dim, div, face, cubes=()):
        """
        dim, div -- as usual
        face  -- pair (coord, value) defining the cube hyperface, with x_{coord}=value
        cubes -- sequence of div-cubes defining a position in the face
        """
        self.dim = dim
        self.div = div
        self.face = tuple(face)
        self.cubes = tuple(tuple(cube) for cube in cubes)

    def is_subset(self, other):
        return self.intersects(other) and len(self.cubes) >= len(other.cubes)

    def intersects(self, other):
        if self.face != other.face:
            return False
        if any(c1 != c2 for c1, c2 in zip(self.cubes, other.cubes)):  # zip shortest
            return False
        return True

    @property
    def face_coord(self):
        return self.face[0]

    @property
    def face_value(self):
        return self.face[1]

    @property
    def _data(self):
        return self.dim, self.div, self.face, self.cubes

    def __eq__(self, other):
        return self._data == other._data

    def __hash__(self):
        return hash(self._data)

    def __lt__(self, other):
        return self._data < other._data

    def map_face(self, cube_map):
        for j, (k, b) in enumerate(cube_map.coords):  # y[j] = x[k]^b
            if k == self.face_coord:
                return j, 1 - self.face_value if b else self.face_value

    def abs_to_face(self, j):
        return j if j < self.face_coord else j-1

    def __rmul__(self, cube_map):
        """Apply cube map (a BaseMap instance)."""
        assert self.dim == cube_map.dim

        new_face = self.map_face(cube_map)

        # j = absolute coordinate, m = relative in self.face, i = relative in new face
        # in-face coordinate map: (u_0,...,u_{face_dim-1}) -> (v_0,...,v_{face_dim-1})
        face_map_data = []
        for j in range(self.dim):
            if j == new_face[0]:
                continue
            k, b = cube_map.coords[j]  # v_i = y_j = x_k^b  (i=index)
            m = self.abs_to_face(k)  # x_k = u_m
            face_map_data.append((m, b))  # v_i = u_m^b

        face_map = BaseMap(face_map_data)
        new_cubes = (face_map.apply_cube(self.div, cube) for cube in self.cubes)
        return type(self)(self.dim, self.div, new_face, new_cubes)

    def all_cube_maps(self):
        if not hasattr(self, '_all_cube_maps'):
            self._all_cube_maps = list(BaseMap.gen_base_maps(self.dim, time_rev=False))
        return self._all_cube_maps

    def argmul_intersect(self, other, bms=None):
        """Yield bms that bm * self intersects other."""
        if bms is None:
            bms = self.all_cube_maps()
        for bm in bms:
            if self.map_face(bm) != other.face:  # optimization!
                continue
            if (bm * self).intersects(other):
                yield bm

    def map_to_cube(self, div, cube):
        assert div == self.div
        assert cube[self.face_coord] == (0 if self.face_value == 0 else div-1)
        face_cube = tuple(cube[c] for c in range(self.dim) if c != self.face_coord)
        return type(self)(self.dim, self.div, self.face, (face_cube,) + self.cubes)

    def gen_neighbours(self):
        if self.face_value == 0:
            delta = -1
            new_face_value = 1
        else:
            delta = 1
            new_face_value = 0
        cube = tuple(0 if j != self.face_coord else delta for j in range(self.dim))
        yield cube, type(self)(self.dim, self.div, (self.face_coord, new_face_value), self.cubes)

    def divide(self, div):
        assert div == self.div
        dim = self.dim
        fc, fv = self.face

        def make_abs(face_cube):
            # face_cube = face_dim-cube on the face; other coords are determined by self.face
            return tuple((0 if fv == 0 else div-1) if j == fc else face_cube[self.abs_to_face(j)] for j in range(dim))

        if self.cubes:
            face_cubes = [self.cubes[0]]
            new_set = type(self)(dim, div, self.face, self.cubes[1:])
        else:
            face_cubes = itertools.product(range(div), repeat=dim-1)
            new_set = self

        for face_cube in face_cubes:
            yield make_abs(face_cube), new_set

    def __str__(self):
        text = '{}={}'.format(self.face_coord, self.face_value)
        if self.cubes:
            cubes_str = '->'.join([str(cube) for cube in self.cubes])
            text += ':{}'.format(cubes_str)
        return text


class Link(namedtuple('Link', ['entrance', 'exit'])):
    """
    Link is defined by pair of subsets of [0,1]^d: entrance and exit subsets.

    Link represents the set of curves with curve.entrance in link.entrance and curve.exit in link.exit.
    Entrance and exit must be subsets of [0,1]^d (derived from Subset)
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

    def __and__(self, other):
        return type(self)(self.entrance & other.entrance, self.exit & other.exit)

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
