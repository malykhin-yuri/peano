"""
This module provides Estimator class to estimate
regular curve or fuzzy curve dilation.
SAT-solvers are used for fuzzy curves.
"""

from collections import Counter, namedtuple
from collections.abc import Sized
from heapq import heappop, heappush
import logging

from quicktions import Fraction

from .utils import get_lcm, get_int_cube_with_cache, get_int_time_with_cache
from . import _sat_adapters
from .curves import Curve


logger = logging.getLogger(__name__)


class _PiecePosition:
    # Time (cnums) and space (cubes) location of a fraction in a curve.
    # A support of the fraction is the last cube in the nested sequence of cubes.

    def __init__(self, dim, div, cnums, cubes):
        self.dim = dim
        self.div = div
        self.cnums = tuple(cnums)
        self.cubes = tuple(cubes)
        self.depth = len(self.cnums)

    def specify(self, cnum, cube):
        return _PiecePosition(
            dim=self.dim,
            div=self.div,
            cnums=self.cnums + (cnum,),
            cubes=self.cubes + (cube,),
        )

    def get_int_coords(self):
        # Natural integer coordinates - time and lower-left cube vertex
        # Returns pair x, t:
        #    x: int cube, cj <= xj <= cj+1; where abs cube: cj/N^l <= xj <= (cj+1)/N^l
        #    t: int time, [t, t+1]; where abs time: [t/G^l, (t+1)/G^l]
        #    These may be viewed as coordinates in curve [0,G^l] -> [0,N^l]^d
        return (
            get_int_cube_with_cache(self.dim, self.div, self.cubes),
            get_int_time_with_cache(self.dim, self.div, self.cnums),
        )


class _CurvePiece(namedtuple('_CurvePiece', ['curve', 'pnum', 'pos'])):
    # Fraction of a curve, defined by a fuzzy curve, pnum and _PiecePosition
    # The specs in the curve must be specified in all cubes of the position, except last one.
    # Helper class, used in _CurvePiecePair only.

    def get_last_spec(self):
        # Get spec X such that the curve on fraction equals X * self.curve.
        # Use this for fully defined fractal curves, because in the fuzzy case
        # the spec on the last cube is not defined.
        return self.curve.get_deep_spec(self.pnum, self.pos.cnums)

    def divide(self):
        # Divide the fraction in all possible ways.
        # define orientation of last_but_one fraction of a curve
        # in the last fraction we do not know the orientation yet!
        prev_spec = self.curve.get_deep_spec(self.pnum, self.pos.cnums[:-1])

        active_pnum = prev_spec.pnum
        active_cnum = self.pos.cnums[-1]  # the cube in prev_curve that will be divided
        orig_cnum = prev_spec.base_map.apply_cnum(self.curve.genus, active_cnum)

        for orig_spec in self.curve.gen_allowed_specs(active_pnum, orig_cnum):
            # we need to get spec in prev_curve to proceed to last curve
            # specs are conjugated, see curves.FuzzyCurve.apply_cube_map
            sp = orig_spec.conjugate_by(prev_spec.base_map)

            specified_curve = self.curve.specify(active_pnum, orig_cnum, orig_spec)

            # last_curve = sp * prev_curve, but we do not use curve mult - optimization
            # last_curve_proto = (sp * prev_curve).proto
            last_curve_proto = (sp.base_map * prev_spec.base_map) * self.curve.patterns[sp.pnum].proto
            for cnum, cube in enumerate(last_curve_proto):
                new_pos = self.pos.specify(cnum, cube)
                new_piece = _CurvePiece(specified_curve, self.pnum, new_pos)
                yield new_piece


class _CurvePiecePair(namedtuple('_CurvePiecePair', ['curve', 'junc', 'pos1', 'pos2'])):
    # Pair or curve fractions.
    # curve: fuzzy curve
    # junc: junction that defines pair
    # pos1: position of first fraction; in the original curve, not rotated(!)
    # pos2: the same for second fraction

    def _get_piece(self, piece_no):
        pnum = self.junc.spec1.pnum if piece_no == 1 else self.junc.spec2.pnum
        pos = self.pos1 if piece_no == 1 else self.pos2
        return _CurvePiece(self.curve, pnum, pos)

    @classmethod
    def init_first_order(cls, curve, junc, cnum1, cnum2):
        # Init a pair or first-order fractions.
        # cnum1, cnum2: cnum of the first/second fraction (in the original curve)
        def get_pos(pnum, cnum):
            cube = curve.patterns[pnum].proto[cnum]
            return _PiecePosition(dim=curve.dim, div=curve.div, cnums=[cnum], cubes=[cube])
        return cls(curve, junc, get_pos(junc.spec1.pnum, cnum1), get_pos(junc.spec2.pnum, cnum2))

    def get_last_specs(self):
        return self._get_piece(1).get_last_spec(), self._get_piece(2).get_last_spec()

    @property
    def depth(self):
        return self.junc.depth + min(self.pos1.depth, self.pos2.depth)

    def divide_balanced(self):
        # Divide one of the fractions keeping pair balanced: depth1 == depth2 or depth1 == depth2 + 1.

        # use curve from divided piece because it has specified curve
        if self.pos1.depth > self.pos2.depth:
            for subpiece in self._get_piece(2).divide():
                yield _CurvePiecePair(subpiece.curve, self.junc, self.pos1, subpiece.pos)
        else:
            for subpiece in self._get_piece(1).divide():
                yield _CurvePiecePair(subpiece.curve, self.junc, subpiece.pos, self.pos2)


class _BoundedItemsHeap:
    # Collection of items with lower/uppper bounds on their "values" (item.lo, item.up)
    #
    # Two thresholds may be set:
    # * good - if item.up <= good then item is considered "good" and discarded
    # * bad - if item.lo >= bad then item is "bad" and it is temporarily stored
    # note that an item may be good and bad simultaneously, if bad_thr <= lo <= up <= good_thr
    # in this case it is considered as "good"
    # * otherwise item is "active", they are stored in heap with priority=up

    def __init__(self, keep_max_lo_item=False, **kwargs):
        # use heapq algorithm for list of tuples (priority, increment, item)
        self._heap = []
        self._bad_items = []
        self._good_threshold = None
        self._bad_threshold = None

        self._keep_max_lo_item = keep_max_lo_item
        if self._keep_max_lo_item:
            self.max_lo_item = None

        self.stats = Counter()
        self.other = kwargs
        self._inc = 0

    def set_good_threshold(self, threshold):
        # If up <= good, item is considered "good" (hence insignificant).
        self.stats['set_good_threshold'] += 1
        self._good_threshold = threshold

    def set_bad_threshold(self, threshold):
        # If lo >= bad, item is considered "bad".
        self.stats['set_bad_threshold'] += 1
        self._bad_threshold = threshold

    def has_items(self):
        return bool(self._heap)

    def items(self):
        # returns iterator over active items
        return (node[-1] for node in self._heap)

    def push(self, item):
        # Add node checking the thresholds:
        # good pair is dropped / bad pair is temporarily stored / otherwise node is added to the heap
        self.stats['push'] += 1

        if self._keep_max_lo_item:
            if self.max_lo_item is None or item.lo > self.max_lo_item.lo:
                self.max_lo_item = item

        # the order of checks may be important
        # is Estimator we add bad pairs to SAT clauses, so it may be more effective to
        # check for goodness first (recall that an item may be good and bad simultaneously)
        if (self._good_threshold is not None) and item.up <= self._good_threshold:
            self.stats['good'] += 1
            return
        if (self._bad_threshold is not None) and item.lo >= self._bad_threshold:
            self._bad_items.append(item)
            self.stats['bad'] += 1
            return

        self._inc += 1
        node = (-item.up, self._inc, item)  # first key is priority; _inc to avoid comparing items
        heappush(self._heap, node)

    def extend(self, items):
        for item in items:
            self.push(item)

    def top(self):
        # Active item with highest priority (up)
        return self._heap[0][-1]

    def pop(self):
        # Pop and return active item with highest priority (up)
        return heappop(self._heap)[-1]

    def pop_bad_items(self):
        # Return bad pairs list and empty it.
        bad_items = self._bad_items
        self._bad_items = []
        return bad_items

    def can_cleanup(self):
        # cleanup pushs <= regular pushs
        return self.stats['cleanup_push'] <= (self.stats['push'] - self.stats['cleanup_push'])

    def cleanup(self):
        # rebuild heap with actual thresholds
        items = list(self.items())
        self._heap = []
        self.extend(items)
        self.stats['cleanup_push'] += len(items)
        self.stats['cleanup_count'] += 1

    def copy_and_cleanup(self, good_threshold=None, bad_threshold=None):
        # copy initial tree and apply thresholds, if any
        assert self._good_threshold is None
        assert self._bad_threshold is None
        assert not self._keep_max_lo_item
        new_heap = _BoundedItemsHeap()
        new_heap.set_good_threshold(good_threshold)
        new_heap.set_bad_threshold(bad_threshold)
        new_heap.extend(self.items())
        new_heap.stats['copy_push'] += len(list(self.items()))
        return new_heap


class Estimator:
    """
    Estimator - estimates curve dilation.

    Dilation of curve gamma: [0,1]->[0,1]^d is defined as
    WD(gamma) := sup_{s,t} ||gamma(s)-gamma(t)||^d / |t-s|

    Attributes:
        stats: counter of some global statistics
    """

    class _RunOutOfIterationsException(Exception):
        pass

    def __init__(self, ratio_func, cache_max_size=2**18, cache_max_depth=4):
        """
        Init Estimator instance.

        Args:
            ratio_func: function (dim, dx, dt) -> Fraction
              it is assumed to be d-uniform and coordinate-monotone
            cache_max_size: cache size for pairs bounds
        """

        self.ratio_func = ratio_func
        self.sum_stats = Counter()
        self.max_stats = {}
        self._get_bounds_cache = {}
        self._cache_max_size = cache_max_size
        self._cache_max_depth = cache_max_depth

    def _get_bounds(self, pair, brkline=None):
        # Get lower and upper bounds for max ratio of given fractions pair:
        #   WD(f1,f2) := sup ||gamma(s)-gamma(t)||^d/|s-t|:  gamma(s) in f1, gamma(t) in f2
        # brkline: instance of _IntegerBrokenLine class (for pcount==1 only!)
        # Returns triple (lo, up, argmax), argmax only for brkline

        self.sum_stats['get_bounds_calls'] += 1
        dim = pair.curve.dim
        N = pair.curve.div

        pos1, pos2 = pair.pos1, pair.pos2
        pair_depth = max(pos1.depth, pos2.depth)  # not count junc.depth
        self.sum_stats['pair_depth'] += pair_depth
        self.max_stats['pair_depth'] = max(self.max_stats.get('pair_depth', 0), pair_depth)

        use_cache = (brkline is None and pair_depth <= self._cache_max_depth)  # not implemented for brkline
        if use_cache:
            # note that pos pnums are not part of the key
            # because cube time/space locations do not depend on them and dilation is not affected
            # however, pnums are important for brkline rotation - but brkline case is not cached
            cache_key = (dim, N, pair.junc, pos1.cnums, pos1.cubes, pos2.cnums, pos2.cubes)
            cache = self._get_bounds_cache
            cache_value = cache.get(cache_key)
            if cache_value is not None:
                self.sum_stats['get_bounds_cache_hit'] += 1
                return cache_value
            else:
                self.sum_stats['get_bounds_cache_miss'] += 1

        # these are integer positions in original curve patterns
        # we will transform them to absolute coords
        x1, t1 = pos1.get_int_coords()
        x2, t2 = pos2.get_int_coords()

        use_brkline = (brkline is not None)
        if use_brkline:
            # rotations for broken lines
            sp1, sp2 = pair.get_last_specs()
            assert sp1.pnum == sp2.pnum == 0
            brk1_bm, brk2_bm = sp1.base_map, sp2.base_map

        junc = pair.junc

        # junc: apply base_maps to coordinates
        pos1_sub_div = N**pos1.depth
        pos1_sub_genus = pos1_sub_div**dim

        pos2_sub_div = N**pos2.depth
        pos2_sub_genus = pos2_sub_div**dim

        jbm1, jbm2 = junc.spec1.base_map, junc.spec2.base_map
        x1 = jbm1.apply_cube(pos1_sub_div, x1)
        t1 = jbm1.apply_cnum(pos1_sub_genus, t1)
        x2 = jbm2.apply_cube(pos2_sub_div, x2)
        t2 = jbm2.apply_cnum(pos2_sub_genus, t2)
        if use_brkline:
            brk1 = ((jbm1 * brk1_bm) * brkline).points
            brk2 = ((jbm2 * brk2_bm) * brkline).points

        # common scale
        if pos1.depth == pos2.depth:
            mx2 = mt2 = 1
        elif pos1.depth == pos2.depth + 1:
            mx2, mt2 = N, N**dim
            x2 = [xj * mx2 for xj in x2]
            t2 *= mt2
            if use_brkline:
                brk2 = [([xj * mx2 for xj in x], t * mt2) for x, t in brk2]
        else:
            raise ValueError("Unbalanced positions!")

        mx = pos1_sub_div
        mt = pos1_sub_genus

        # now we have the following integer coordinates:
        # cube1: x1j <= xj <= x1j + 1    -- cube inside [0, mx]^d
        # cube2: x2j <= xj <= x2j + mx2  -- cube inside [0, mx2 * mx]^d, + shift junc_dx * mx
        #
        # time1: t1 <= t <= t1 + 1
        # time2: t2 <= t <= t2 + mt2,  + shift junc_dt * mt

        # junc: shifts
        t2 += junc.delta_t * mt
        x2 = [x2j + dxj * mx for x2j, dxj in zip(x2, junc.delta_x)]

        max_dx = [max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2)) for x1j, x2j in zip(x1, x2)]

        max_dt = t2 + mt2 - t1  # max(t_2 - t_1)
        min_dt = t2 - (t1 + 1)  # min(t_2 - t_1)

        # also ratio(min_dx, min_dt) is a lower bound, but inefficient
        lo = self.ratio_func(dim, max_dx, max_dt)
        up = self.ratio_func(dim, max_dx, min_dt)

        argmax = None
        if use_brkline:
            brk_mx, brk_mt = brkline.mx, brkline.mt
            x1 = [xj * brk_mx for xj in x1]
            x2 = [xj * brk_mx for xj in x2]
            t1 *= brk_mt
            t2 *= brk_mt
            for x1rel, t1rel in brk1:
                x1_point = [x1j + x1relj for x1j, x1relj in zip(x1, x1rel)]
                t1_point = t1 + t1rel
                for x2rel, t2rel in brk2:
                    x2_point = [x2j + x2relj for x2j, x2relj in zip(x2, x2rel)]
                    t2_point = t2 + t2rel

                    dx = [x1j - x2j for x1j, x2j in zip(x1_point, x2_point)]
                    dt = t2_point - t1_point
                    lo_point = self.ratio_func(dim, dx, dt)

                    if lo_point > lo or argmax is None:
                        lo = lo_point
                        x1_real = [Fraction(x1j, mx * brk_mx) for x1j in x1_point]
                        x2_real = [Fraction(x2j, mx * brk_mx) for x2j in x2_point]
                        t1_real = Fraction(t1_point, mt * brk_mt)
                        t2_real = Fraction(t2_point, mt * brk_mt)
                        argmax = {'x1': x1_real, 't1': t1_real, 'x2': x2_real, 't2': t2_real, 'junc': junc}

        result = (lo, up, argmax)
        if use_cache:
            self.max_stats['cache_size'] = max(self.max_stats.get('cache_size', 0), len(cache))
            if self._cache_max_size is not None and len(cache) >= self._cache_max_size:
                # poor man's LRU cache :(
                cache.clear()
                self.sum_stats['get_bounds_cache_cleanup'] += 1
            cache[cache_key] = result

        return result

    def _create_tree(self, curve, good_threshold=None, bad_threshold=None, **tree_kwargs):
        # Create initial pairs tree from a curve
        G = curve.genus
        tree = _BoundedItemsHeap(**tree_kwargs)
        tree.set_good_threshold(good_threshold)
        tree.set_bad_threshold(bad_threshold)
        for junc in curve.gen_auto_junctions():
            for cnum1 in range(G):
                for cnum2 in range(cnum1 + 2, G):
                    pair = _CurvePiecePair.init_first_order(curve, junc, cnum1, cnum2)
                    self._push_tree(tree, pair)

        for junc in curve.gen_regular_junctions():
            last_cnum1 = 0 if junc.spec1.base_map.time_rev else G - 1
            first_cnum2 = G - 1 if junc.spec2.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) == (last_cnum1, first_cnum2):
                        continue
                    pair = _CurvePiecePair.init_first_order(curve, junc, cnum1, cnum2)
                    self._push_tree(tree, pair)

        return tree

    # items for heap = pairs of curve fractions ( _CurvePiecePair)
    # together with lo/up bounds on their dilation
    _BoundedPair = namedtuple('_BoundedPair', ['pair', 'lo', 'up', 'argmax'])

    def _push_tree(self, tree, pair):
        lo, up, argmax = self._get_bounds(pair, brkline=tree.other.get('brkline'))
        item = self._BoundedPair(pair, lo, up, argmax)
        tree.push(item)

    def _divide_tree(self, tree):
        # Divide worst pair
        worst_item = tree.pop()
        for new_pair in worst_item.pair.divide_balanced():
            self._push_tree(tree, new_pair)

    def _add_tree_stats(self, tree):
        self.sum_stats.update({'ptree.{}'.format(k): v for k, v in tree.stats.items()})

    def _update_max_stats(self, k, v):
        max_stats = self.max_stats
        if k in max_stats:
            max_stats[k] = max(max_stats[k], v)
        else:
            max_stats[k] = v

    def _forbid(self, tree, adapter):
        # forbid bad configurations and drop them from the tree
        for item in tree.pop_bad_items():
            adapter.add_forbid_clause(item.pair.junc, item.pair.curve)

    @staticmethod
    def _try_by_strategy(info):
        strategy = info['strategy']
        info.setdefault('iter', 0)
        info['iter'] += 1
        if strategy['type'] == 'equal':
            if info['iter'] % strategy['count'] == 0:
                return True
        elif strategy['type'] == 'geometric':
            info.setdefault('next_try_iter', 1)
            if info['iter'] >= info['next_try_iter']:
                info['next_try_iter'] = int(info['next_try_iter'] * strategy['multiplier']) + 1
                return True

    def estimate_dilation(self, curve, *args, **kwargs):
        """Dispatcher method: uses estimate_dilation_regular or estimate_dilation_fuzzy"""
        if isinstance(curve, Curve):
            return self.estimate_dilation_regular(curve, *args, **kwargs)
        else:
            return self.estimate_dilation_fuzzy(curve, *args, **kwargs)

    def estimate_dilation_regular(self, curve, rel_tol_inv=100, max_iter=None, use_vertex_brkline=False, max_depth=None):
        """
        Estimate dilation for a regular peano curve (class Curve).

        Args:
            curve: Curve instance, fully defined polyfractal curve
            rel_tol_inv: inverted relative tolerance (may be set to None)
            max_iter: limit for subdivisions iters
            use_vertex_brkline: use vertex moments (broken line) for dilation lower bounds
            max_depth: allow not to consider pairs of higher depth
              note that fraction depth = junc depth + piece depth;
              if dilation is attained at pair of fractions of depth <= max_depth,
              then resulting "lo" will be exact

        Returns:
            dict with keys:
            'lo': lower bound - dilation(curve) >= lo
            'up': upper bound - dilation(curve) <= up
            'argmax': pair of points where lo is achieved (if use_vertex_brkline is set)
        """

        tree_kwargs = {'keep_max_lo_item': True}
        if use_vertex_brkline:
            if curve.pcount > 1:
                raise NotImplementedError("Brklines for multiple patterns not implemented!")
            vertex_brkline = list(curve.get_vertex_moments().items())
            tree_kwargs['brkline'] = _IntegerBrokenLine.init_from_rational(curve.dim, vertex_brkline)

        if rel_tol_inv is not None:
            tolerance = Fraction(rel_tol_inv + 1, rel_tol_inv)

        # This is a basic algorithm that does not require SAT-solvers.
        # We maintain a "tree" (_BoundedItemsHeap instance) of pairs (_BoundedPair)
        # of all non-adjacent curve fractions, for all junctions.
        # For each pair we get lower and upper bounds for dilation (see _get_bound)
        # At each iteration we divide the worst pair (with max upper bound)
        # Dilation is attended at one of the active pairs of the tree.

        pairs_tree = self._create_tree(curve, **tree_kwargs)

        # if there is a pair with dilation >= lo, then lo is the bound for the whole curve
        # and we do not need to consider pairs with less-or-equal upper bound
        curr_lo = pairs_tree.max_lo_item.lo
        argmax = pairs_tree.max_lo_item.argmax
        pairs_tree.set_good_threshold(curr_lo)

        # upper bound for worst active pair gives us upper bound for the curve
        curr_up = pairs_tree.top().up
        logger.info('start bounds: %.5f < %.5f', curr_lo, curr_up)

        # invariant: dilation of the curve is in [curr_lo, curr_up]
        iter_no = 0
        while True:
            iter_no += 1
            if (rel_tol_inv is not None) and curr_up <= curr_lo * tolerance:
                break
            if (max_iter is not None) and iter_no > max_iter:
                break
            if max_depth is not None:
                # given max_depth, we proceed as usual and consider deep fractions also
                # to obtain good up bounds and truncate tree faster;
                # we wait till all active pairs become deep enough to ensure that
                # all pairs with depth <= max_depth had been processed;
                # but we have to cleanup to discard old pairs with not-so-high up
                if pairs_tree.can_cleanup():
                    pairs_tree.cleanup()
                    depth = min(item.pair.depth for item in pairs_tree.items())
                    if depth > max_depth:
                        break

            self._divide_tree(pairs_tree)

            item = pairs_tree.max_lo_item
            if item.lo > curr_lo:
                curr_lo = item.lo
                argmax = item.argmax
                pairs_tree.set_good_threshold(curr_lo)
                logger.info('new lower bound: %.5f < %.5f', curr_lo, curr_up)

            new_up = pairs_tree.top().up
            if new_up < curr_up:
                logger.info('new upper bound: %.5f < %.5f', curr_lo, curr_up)
                curr_up = new_up

        self._add_tree_stats(pairs_tree)
        res = {'up': curr_up, 'lo': curr_lo}
        if argmax is not None:
            res['argmax'] = argmax

        return res

    def bisect_dilation_fuzzy(self, curve, good_threshold, bad_threshold,
                              max_iter=None, sat_strategy=None,
                              init_pairs_tree=None, init_sat_adapter=None):
        """
        Decide if there is a "good" regular curve or all curves are "bad".

        Consider the minimal dilation WD(gamma) for all regular curves from given fuzzy curve.
        Given two thresholds: bad_thr < good_thr, method either:
        * finds a "good" (dilation <= good_thr) regular curve from this fuzzy curve; returns that curve
          so, min WD <= good_thr
        * proves that dilation is "bad" (>= bad_thr) for all curves; returns no curve
          so, min WD >= bad_thr in this case
        If (bad_thr <= min WD <= good_thr), both cases take place and return value is not specified,
        so the resulting curve is guaranteed only if min WD < bad_thr.

        Args:
            curve: FuzzyCurve instance
            good_threshold: good curves are those with dilation <= good_thr
            bad_threshold: bad curves are those with dilation >= bad_thr
              It is required that bad_threshold < good_threshold
            max_iter: max number of divisions; raise exception if maximum iterations reached
            sat_strategy: when do we call sat solver:
              strategy['type'] == 'equal': call every strategy['count'] divisions (default)
              strategy['type'] == 'geometric': call on strategy['multiplier']**k iterations
            init_pairs_tree: cached _create_tree
            init_sat_adapter: cached initial sat_adapter

        Returns:
            curve or None (if not found)
        """

        # This is the main algorithm of the whole "peano" package.
        #
        # All possible regular curves from given fuzzy curve are encoded using
        # boolean variables in sat adapter.
        #
        # We grow the pairs tree with fixed good and bad thresholds
        # all bad pairs are added to the list of forbidden configurations,
        # i.e. a list of boolean clauses in sat adapter.
        #
        # If we can't find a model, then all curves are necessarily bad
        # (i.e., all possible configurations are forbidden)
        # If there are no active pairs in the tree and adapter finds a model,
        # then the corresponding curve is good because it avoids bad pairs
        # and all other curve's pairs are good


        # how often should we call sat solver? default is equidistant strategy
        if sat_strategy is None:
            sat_strategy = {'type': 'equal', 'count': 100}
        try_sat_info = {'strategy': sat_strategy}

        thrs = {'good_threshold': good_threshold, 'bad_threshold': bad_threshold}
        if init_pairs_tree is None:
            pairs_tree = self._create_tree(curve, **thrs)
        else:
            pairs_tree = init_pairs_tree.copy_and_cleanup(**thrs)

        if init_sat_adapter is None:
            adapter = _sat_adapters.CurveSATAdapter(curve)
        else:
            adapter = init_sat_adapter.copy()

        self._forbid(pairs_tree, adapter)

        sum_stats = Counter()
        no_model = None
        iter_no = 0
        while pairs_tree.has_items():
            iter_no += 1
            if (max_iter is not None) and iter_no > max_iter:
                raise self._RunOutOfIterationsException()

            self._divide_tree(pairs_tree)
            self._forbid(pairs_tree, adapter)
            self.sum_stats['divide_iter'] += 1

            try_sat = self._try_by_strategy(try_sat_info)
            if not try_sat:
                continue

            logger.debug('iter %d, try SAT solver; stats: %s, %s', iter_no, self.sum_stats, self.max_stats)
            self.sum_stats['sat.solve_calls'] += 1
            if not adapter.solve():
                no_model = True
                break

        self._add_tree_stats(pairs_tree)
        self._update_max_stats('sat.problem_size', adapter.get_problem_size())

        if no_model or not adapter.solve():
            return None

        model = adapter.get_model()
        return adapter.get_curve_from_model(model)

    def estimate_dilation_fuzzy(self, curve, rel_tol_inv=1000, stop_upper_bound=None,
                                start_lower_bound=None, start_upper_bound=None, start_curve=None,
                                **kwargs):
        """
        Estimate minimal dilation of a fuzzy curve.

        We estimate min WD(gamma) for regular curves gamma from given fuzzy curve,
        i.e., the dilation of best curve.

        Args:
            curve: fuzzy curve
            rel_tol_inv: inverted relative tolerance, integer
            stop_upper_bound: do not proceed if min WD is higher than this

            start_lower_bound: known lower bound on min WD, start bisection with it
            start_upper_bound: known upper bound on min WD, start bisection with it
            start_curve: known curve with start_upper_bound

            **kwargs: passed to bisect_dilation_fuzzy

        Returns:
             dict with keys:
            'lo': lower_bound
            'up': upper_bound
            'curve': curve_example with dilation in [lo, up]
        """

        # This method is simply "bisection" algorithm based on bisect_dilation_fuzzy.

        if start_lower_bound is None:
            curr_lo = Fraction(0)
        else:
            curr_lo = start_lower_bound

        if start_upper_bound is None:
            # got from first regular curve
            curve0 = curve.get_curve_example()
            curr_up = self.estimate_dilation_regular(curve0, rel_tol_inv=rel_tol_inv)['up']
            curr_curve = curve0
        else:
            curr_up = start_upper_bound
            curr_curve = start_curve

        init_pairs_tree = self._create_tree(curve)
        self._add_tree_stats(init_pairs_tree)  # will not use init_pairs_tree anymore

        init_sat_adapter = _sat_adapters.CurveSATAdapter(curve)

        # invariants:
        # * minimum dilation is in [curr_lo, curr_up]
        # * curr_curve dilation also in [curr_lo, curr_up]
        tolerance = Fraction(rel_tol_inv + 1, rel_tol_inv)
        while curr_up > curr_lo * tolerance:
            self.sum_stats['bisect_iter'] += 1

            if curr_lo == Fraction(0):
                # optimize, 0 is too rude lower bound
                test_lo = Fraction(1, 2) * curr_up
                test_up = Fraction(2, 3) * curr_up
            else:
                test_lo = Fraction(2, 3) * curr_lo + Fraction(1, 3) * curr_up
                test_up = Fraction(1, 3) * curr_lo + Fraction(2, 3) * curr_up

            logger.info(
                'Bisect #%d. best in: [%.5f, %.5f]; seek with thresholds: [%.5f, %.5f]', self.sum_stats['bisect_iter'],
                curr_lo, curr_up, test_lo, test_up,
            )
            logger.debug('precise test thresholds: %s, %s', test_lo, test_up)
            try:
                bisect_result = self.bisect_dilation_fuzzy(
                    curve,
                    bad_threshold=test_lo,
                    good_threshold=test_up,
                    init_pairs_tree=init_pairs_tree,
                    init_sat_adapter=init_sat_adapter,
                    **kwargs,
                )
            except self._RunOutOfIterationsException:
                # run out of iterations
                break

            if bisect_result is not None:
                curr_curve = bisect_result
                curr_up = test_up
            else:
                curr_lo = test_lo

            if (stop_upper_bound is not None) and curr_lo > stop_upper_bound:
                break

        return {
            'curve': curr_curve,
            'lo': curr_lo,
            'up': curr_up,
        }

    def estimate_dilation_sequence(self, curves, rel_tol_inv=1000, rel_tol_inv_mult=4, upper_bound=None, **kwargs):
        """
        Estimate minimal dilation for sequence of fuzzy curves.

        We estimate min_{fuzzy} min_{curve in fuzzy} WD(curve).

        Args:
            curves: iterable of fuzzy curves
            upper_bound: apriori upper bound for best ratio (may return None if violated)
            rel_tol_inv: target tolerance
            rel_tol_inv_mult: current rel_tol_inv is multiplied by this every epoch
            **kwargs: passed as is to estimate_dilation_fuzzy

        Returns:
            dict with keys:
            'lo', 'up' - bounds for minimal dilation
            'curves' - list of regular curves in [lo, up]
            'idxs' - list of corresponding indices of input
        """

        # This method is based on estimate_dilation_fuzzy.
        # We iterate over curves many times with increasing up/lo estimation tolerance.

        tolerance = Fraction(rel_tol_inv + 1, rel_tol_inv)
        CurveItem = namedtuple('CurveItem', 'priority idx lo up curve example'.split())

        def get_item(curve, idx, lo=None, up=None, example=None):
            # store idx in second field to avoid lo/up/... comparison
            priority = -lo if lo is not None else None
            return CurveItem(priority, idx, lo, up, curve, example)

        curr_lo = Fraction(0)
        curr_up = upper_bound

        # Invariants:
        # * minimal dilation is in [curr_lo, curr_up]
        # (curr_up will be defined after first call to estimate_dilation_fuzzy)
        # * curve with min dilation is in active list

        active = (get_item(curve, idx) for idx, curve in enumerate(curves))
        if isinstance(curves, Sized):
            active = list(active)
        curr_rel_tol_inv = 1
        epoch = 0
        while curr_up is None or curr_up > curr_lo * tolerance:
            curr_rel_tol_inv *= rel_tol_inv_mult
            epoch += 1
            total = len(active) if isinstance(active, list) else -1
            new_active = []  # heap of CurveItem
            for cnt, item in enumerate(active):
                if epoch == 1:
                    self.sum_stats['seen_pcurve'] += 1
                logger.info('Epoch #%d, curve %d / %d', epoch, cnt + 1, total)
                res = self.estimate_dilation_fuzzy(
                    item.curve, rel_tol_inv=curr_rel_tol_inv, stop_upper_bound=curr_up,
                    start_lower_bound=item.lo, start_upper_bound=item.up, start_curve=item.example,
                    **kwargs,
                )
                if curr_up is None or res['up'] < curr_up:
                    curr_up = res['up']
                    logger.info('new upper bound: %.3f', curr_up)

                if res['lo'] <= curr_up:
                    # has a chance to be the best
                    new_item = get_item(
                        idx=item.idx,
                        curve=item.curve,
                        lo=res['lo'], up=res['up'],
                        example=res['curve'],
                    )
                    heappush(new_active, new_item)
                    logger.info('added new active item!')

                # it is important to cleanup new_active (in first run) to reduce memory consumption
                while new_active and new_active[0].lo > curr_up:  # priority = -lo
                    heappop(new_active)

                logger.info('current active: %d', len(new_active))

            if not new_active:
                # upper bound is too strong
                return None

            active = sorted(new_active, key=lambda item: item.up)  # better to start with good curves
            curr_lo = min(item.lo for item in active)
            logger.info('current bounds: [%.5f, %.5f]', curr_lo, curr_up)

        return {
            'lo': curr_lo, 'up': curr_up,
            'curves': [item.example for item in active],
            'idxs': [item.idx for item in active],
        }


class _IntegerBrokenLine(namedtuple('_IntegerBrokenLine', ['mx', 'mt', 'points'])):
    @classmethod
    def init_from_rational(cls, dim, brkline):
        denoms = set()
        for x, t in brkline:
            if isinstance(t, Fraction):
                denoms.add(t.q)
            for xj in x:
                if isinstance(xj, Fraction):
                    denoms.add(xj.q)
        lcm = get_lcm(denoms)
        mx = lcm
        mt = lcm**dim
        points = []
        for x, t in brkline:
            xp = tuple(int(Fraction(xj) * mx) for xj in x)
            tp = int(Fraction(t) * mt)
            points.append((xp, tp))
        return cls(mx, mt, points)

    def __rmul__(self, base_map):
        points = [(base_map.apply_x(x, mx=self.mx), base_map.apply_t(t, mt=self.mt)) for x, t in self.points]
        return type(self)(self.mx, self.mt, points)
