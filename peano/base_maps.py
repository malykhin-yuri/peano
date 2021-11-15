import itertools

from quicktions import Fraction


class BaseMap:
    """
    Base map: isometry of cube and (possibly) time reversal.

    Immutable and hashable.
    Acts on function f:[0,1]->[0,1]^d as  Bf: [0,1]--B_time-->[0,1]--f-->[0,1]^d--B_cube-->[0,1]^d

    We severily use caching, because there are not so many possible base maps.

    TODO: maybe, totally cache it for dim <= 4 ?
    """

    basis_letters = 'ijklmn'

    _obj_cache = {}

    def __new__(cls, coords, time_rev=False):
        """
        Create a BaseMap instance.

        Cached method. 
        Params:
            coords:     list of pairs (k, b) defining the isometry
                        B_cube(x_0,...,x_{d-1}) = (y_0,...,y_{d-1}),
                        where y_j = 1 - x_{k_j} if b_j else x_{k_j}

            time_rev:   time reversal (boolean), default: False
                        B_time(t) = 1-t if time_rev else t
        """
        
        coords = tuple((k, bool(b)) for k, b in coords)
        time_rev = bool(time_rev)

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
            coords = [(k, False) for k in range(dim)]
            bm = cls._id_map_cache[dim] = cls(coords)
        return bm

    @classmethod
    def parse(cls, text):
        """
        Convenient way to represent a base map.

        basis -- string, e.g., 'Ij', representing images of standard basis e_1, e_2, ...
            i=e_1, j=e_2, k=e_3, l=e_4, m=e_5, n=e_6
            upper-case letters correspond to negation of vectors: I = -i, J = -j, ...
            To get time reverse, place '~' at the end of the string, e.g., 'iJ~'

        We restrict here to dim <= 6, but the limit may by increased by extending basis_letters.
        """

        text = text.strip()
        if text[-1] == '~':
            time_rev = True
            basis = text[:-1]
        else:
            time_rev = False
            basis = text

        assert len(basis) <= len(cls.basis_letters)

        l2i = {l: i for i, l in enumerate(cls.basis_letters)}
        cmap = {}
        for k, l in enumerate(basis):
            ll = l.lower()
            i, b = l2i[ll], (l != ll)
            # e_k -> (+/-)e_i, so y_i = x_k^b
            cmap[i] = (k, b)

        coords = [cmap[i] for i in range(len(cmap))]
        return cls(coords, time_rev)

    def cube_map(self):
        """Isometry without time reversal."""
        return BaseMap(self.coords)

    def __eq__(self, other):
        return self._key == other._key

    def __hash__(self):
        return self._hash

    def __str__(self):
        bases = [None] * self.dim
        for i, (k, b) in enumerate(self.coords):  # y_i = x_k^b
            img = self.basis_letters[i]
            bases[k] = img.upper() if b else img
        return ''.join(bases) + ('~' if self.time_rev else '')

    def __mul__(self, other):
        """
        Composition of base maps: A * B.
        """
        if not isinstance(other, BaseMap):
            return NotImplemented

        key = other._key
        val = self._mul_cache.get(key, None)
        if val is None:
            # actual multiplication
            assert self.dim == other.dim
            coords = []
            for i in range(self.dim):
                p, b1 = self.coords[i]
                k, b2 = other.coords[p]
                coords.append((k, b1 ^ b2))
            time_rev = self.time_rev ^ other.time_rev
            val = BaseMap(coords, time_rev)
            self._mul_cache[key] = val
        return val

    def get_inverse(self):
        """Inverse of base map B: such X that B*X=X*B=id."""
        val = self._inv_cache
        if val is None:
            # actual inversion
            coords = [None] * self.dim
            for i, (k, b) in enumerate(self.coords):
                coords[k] = (i, b)
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
        """Conjugation."""
        return other * self * other.get_inverse()

    def apply_x(self, x, mx=1):
        """Apply isometry to a point x of [0,1]^d."""
        return tuple(mx-x[k] if b else x[k] for k, b in self.coords)

    def apply_t(self, t, mt=1):
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
        time_rev_variants = [True, False] if time_rev is None else [time_rev]
        for perm in itertools.permutations(range(dim)):
            for flip in itertools.product([True, False], repeat=dim):
                for time_rev in time_rev_variants:
                    yield cls(zip(perm, flip), time_rev)

    @classmethod
    def gen_constraint_fast(cls, src, dst):
        """Generate bms: bm*src == dst. Used for rationals."""
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

        bvars = []
        for xj in src:
            bvar = []
            if xj in group_id:
                bvar.append((0, group_id[xj]))
            xj_compl = Fraction(1) - xj
            if xj_compl in group_id:
                bvar.append((1, group_id[xj_compl]))
            bvars.append(bvar)

        # TODO: try to use combinations_product
        for bvar_list in itertools.product(*bvars):
            good = True
            gcnt = {}
            bs = [h[0] for h in bvar_list]  # list of b's
            gs = [h[1] for h in bvar_list]  # list of group ids
            ps = []  # list of places in group
            for g in gs:
                if g not in gcnt:
                    gcnt[g] = 0
                else:
                    gcnt[g] += 1
                if gcnt[g] >= len(group_coords[g]):
                    good = False
                    break  # too many in group
                ps.append(gcnt[g])

            if not good:
                continue

            for perms in group_perm_list:
                coords = [None] * dim
                for k in range(dim):
                    new_k = perms[gs[k]][ps[k]]  # x_k goes to x_{new_k}
                    coords[new_k] = (k, bs[k])
                yield cls(coords)


    @classmethod
    def std(cls, x):
        """Put x in "standard" (minimal) position, return (std_x, std_bm)."""
        dim = len(x)
        minx, minbm = None, None
        for bm in cls.gen_base_maps(dim, time_rev=False):
            bmx = bm * x
            if minx is None or bmx < minx:
                minx = bmx
                minbm = bm
        return minx, minbm


class Spec:
    """
    Specification of poly-fractal curve to a fraction: BaseMap and pattern choice.

    Specs act on polyfractal curves and form a semi-group.
    """

    # we do not make obj_id, to provide consistency for pickling
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
