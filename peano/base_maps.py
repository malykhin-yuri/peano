import itertools
from collections import Counter

from quicktions import Fraction

from .utils import BASIS_LETTERS


class BaseMap:
    """
    Base map: isometry of cube and (possibly) time reversal.

    Immutable and hashable.
    Acts on function f:[0,1]->[0,1]^d as  Bf: [0,1]--B_time-->[0,1]--f-->[0,1]^d--B_cube-->[0,1]^d
    """

    _obj_cache = {}

    def __new__(cls, coords, time_rev=False):
        """
        Create a BaseMap instance.
        Args:
            coords: tuple of pairs (k, b) defining the isometry
              B_cube(x_0,...,x_{d-1}) = (y_0,...,y_{d-1}),
              where y_j = 1 - x_{k_j} if b_j else x_{k_j}
            time_rev: time reversal (boolean), default: False
              B_time(t) = 1-t if time_rev else t
        """
        coords = tuple(coords)
        time_rev = bool(time_rev)

        # use total object caching, because there are not so many possible base maps.
        # do not make obj_id to provide consistency for pickling
        cache = cls._obj_cache
        key = (coords, time_rev)
        if key in cache:
            return cache[key]

        obj = super().__new__(cls)
        obj.dim = len(coords)
        obj.coords = coords
        obj.time_rev = time_rev
        obj._key = key
        obj._mul_cache = {}
        obj._inv_cache = None
        obj._hash = hash(key)

        cache[key] = obj
        return obj

    _id_map_cache = {}
    @classmethod
    def id_map(cls, dim):
        """Identity map."""
        bm = cls._id_map_cache.get(dim)
        if bm is None:
            coords = tuple((k, False) for k in range(dim))
            bm = cls._id_map_cache[dim] = cls(coords)
        return bm

    @classmethod
    def parse(cls, text):
        """
        Get base map from a convenient string representation.

        Args:
            text: string, representing images of standard basis e_1, e_2, ...
              i=e_1, j=e_2, k=e_3, l=e_4, m=e_5, n=e_6
              upper-case letters correspond to negation of vectors: I = -i, J = -j, ...
              To get time reverse, place '~' at the end of the string, e.g., 'iJ~'
              e.g., 'Ji~' means that e_1->-e_2, e_2->e_1, time rev; (x,y)->(y,-x),t->1-t

        We restrict here to dim <= 6, but the limit may by increased by extending basis_letters.
        """
        text = text.strip()
        if text[-1] == '~':
            time_rev = True
            basis = text[:-1]
        else:
            time_rev = False
            basis = text

        assert len(basis) <= len(BASIS_LETTERS)

        l2i = {l: i for i, l in enumerate(BASIS_LETTERS)}
        cmap = {}
        for k, l in enumerate(basis):
            i = l2i[l.lower()]
            b = (l != l.lower())
            # e_k -> (+/-)e_i, so y_i = x_k^b
            cmap[i] = (k, b)

        coords = tuple(cmap[i] for i in range(len(cmap)))
        return cls(coords, time_rev)

    def cube_map(self):
        """Isometry without time reversal."""
        return BaseMap(self.coords)

    def __eq__(self, other):
        return self._key == other._key

    def __hash__(self):
        return self._hash

    def __str__(self):
        basis = [None] * self.dim
        for i, (k, b) in enumerate(self.coords):  # y_i = x_k^b
            img = BASIS_LETTERS[i]
            basis[k] = img.upper() if b else img
        return ''.join(basis) + ('~' if self.time_rev else '')

    def __mul__(self, other):
        """Composition of base maps: A * B."""
        if not isinstance(other, BaseMap):
            return NotImplemented

        key = other._key
        val = self._mul_cache.get(key)
        if val is None:
            assert self.dim == other.dim
            coords = [None] * self.dim
            for i in range(self.dim):
                p, b1 = self.coords[i]  # A: y->z; z_i = y_p^b1
                k, b2 = other.coords[p]  # B: x->y; y_p = x_k^b2
                coords[i] = (k, b1 ^ b2)  # AB: x->y; z_i = x_k^(b1 xor b2)
            self._mul_cache[key] = val = BaseMap(coords, self.time_rev ^ other.time_rev)
        return val

    def get_inverse(self):
        """Group inverse of base map B: such X that B*X=X*B=id."""
        val = self._inv_cache
        if val is None:
            coords = [None] * self.dim
            for i, (k, b) in enumerate(self.coords):
                coords[k] = (i, b)  # y_i = x_k^b <=> x_k = y_i^b
            self._inv_cache = val = BaseMap(coords, self.time_rev)
        return val

    def __pow__(self, power):
        if power < 0:
            return self.get_inverse()**(-power)
        elif power == 1:
            return self
        elif power % 2 == 0:
            p1 = power // 2
            t = self**p1
            return t*t
        else:
            p1 = (power - 1) // 2
            t = self**p1
            return t*t*self

    def __invert__(self):
        """Time reversion, do not confuse with group inverse."""
        return BaseMap(self.coords, not self.time_rev)

    def conjugate_by(self, other):
        """Conjugation: g * X * g^{-1}."""
        return other * self * other.get_inverse()

    def apply_x(self, x, mx=Fraction(1)):
        """Apply cube isometry to a point x of [0,mx]^d."""
        return tuple(mx-x[k] if b else x[k] for k, b in self.coords)

    def apply_t(self, t, mt=Fraction(1)):
        """Apply time isometry to a point t of [0,mt]."""
        return mt - t if self.time_rev else t

    def apply_vec(self, v):
        """Apply linear part of isometry to a vector."""
        return tuple(-v[k] if b else v[k] for k, b in self.coords)

    def apply_cube(self, div, cube):
        """Apply isometry to a sub-cube."""
        return tuple(div-cube[k]-1 if b else cube[k] for k, b in self.coords)

    def apply_cube_start(self, cube_start, cube_length):
        """Apply isometry to a cube of given length, return it's min (start) vertex."""
        return tuple(1-cube_start[k]-cube_length if b else cube_start[k] for k, b in self.coords)

    def apply_cnum(self, genus, cnum):
        return genus - 1 - cnum if self.time_rev else cnum

    @classmethod
    def gen_base_maps(cls, dim, time_rev=None):
        """
        Generate all base maps of given dimension.

        Args:
            time_rev: if not None, yield only bms with given time_rev
        """
        time_rev_variants = (True, False) if time_rev is None else (time_rev,)
        for perm in itertools.permutations(range(dim)):
            for flip in itertools.product([True, False], repeat=dim):
                for time_rev in time_rev_variants:
                    yield cls(zip(perm, flip), time_rev)

    @classmethod
    def gen_constraint_fast(cls, src, dst):
        """
        Generate cube maps that map src to dst.

        Args:
            src, dst: points, i.e. iterables of fractions (dim=length)

        Yields:
            bm: bm*src == dst, bm without time_rev
        """
        dim = len(src)
        group_id = {}
        group_coords = []
        for k, xj in enumerate(dst):
            gid = group_id.get(xj)
            if gid is None:
                gid = len(group_id)
                group_id[xj] = gid
                group_coords.append([])
            group_coords[gid].append(k)

        group_perms = [itertools.permutations(gc) for gc in group_coords]
        group_perm_list = list(itertools.product(*group_perms))

        flip_variants = []  # possible variants (b, group_id) for each coordinate
        for xj in src:
            var = []
            if xj in group_id:
                var.append((False, group_id[xj]))
            xj_compl = Fraction(1) - xj
            if xj_compl in group_id:
                var.append((True, group_id[xj_compl]))
            flip_variants.append(var)

        for flip_list in itertools.product(*flip_variants):
            good = True
            gcnt = Counter()

            # try to check if there exists a bm that maps src->dst
            for _, g in flip_list:
                gcnt[g] += 1
                if gcnt[g] > len(group_coords[g]):
                    good = False
                    break  # too many coords get in group, no bms for this flip_list

            if not good:
                continue

            for perms in group_perm_list:
                perms = tuple(list(perm) for perm in perms)
                coords = [None] * dim
                for k, (b, g) in enumerate(flip_list):
                    new_k = perms[g].pop(0)
                    coords[new_k] = (k, b)
                yield cls(coords)


class Spec:
    """
    Specification of poly-fractal curve to a fraction: BaseMap and pattern choice.

    Specs act on polyfractal curves and form a semi-group,
    i.e. Sf=(bm,pnum)f means bm-mapped curve f with pnum pattern selected.
    """

    _obj_cache = {}

    def __new__(cls, base_map, pnum=0):
        key = (base_map._key, pnum)
        cache = cls._obj_cache
        if key in cache:
            return cache[key]

        obj = super().__new__(cls)
        obj.base_map = base_map
        obj.pnum = pnum
        obj._key = key
        obj._hash = hash(key)
        cache[key] = obj
        return obj

    def __eq__(self, other):
        return self._key == other._key

    def __hash__(self):
        return self._hash

    @classmethod
    def parse(cls, text):
        """Parse spec, e.g. 3ijk for pnum=3 and bm=ijk."""
        if text[0].isdigit():
            pnum = int(text[0])
            basis = text[1:]
        else:
            pnum = 0
            basis = text
        return cls(base_map=BaseMap.parse(basis), pnum=pnum)

    def __str__(self):
        return '{}{}'.format(self.pnum, self.base_map)

    def __mul__(self, other):
        """
        Composition of specs, i.e., (S1*S2) f = S1(S2 f).

        We also allow to combine specs and base maps,
        as they act on curves too.
        """
        if isinstance(other, BaseMap):
            other_bm = other
        elif isinstance(other, Spec):
            other_bm = other.base_map
        else:
            return NotImplemented
        return Spec(self.base_map * other_bm, self.pnum)

    def __rmul__(self, other):
        if isinstance(other, BaseMap):
            other_bm = other
            pnum = self.pnum
        elif isinstance(other, Spec):
            other_bm = other.base_map
            pnum = other.pnum
        else:
            return NotImplemented
        return Spec(other_bm * self.base_map, pnum)

    def __invert__(self):
        return Spec(~self.base_map, self.pnum)

    def conjugate_by(self, other):
        """Conjugate by base map (not spec!)."""
        assert isinstance(other, BaseMap)
        return Spec(self.base_map.conjugate_by(other), self.pnum)
