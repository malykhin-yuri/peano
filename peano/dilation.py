"""
This module provides Estimator class to estimate
regular curve or fuzzy curve dilation.
SAT-solvers are used for fuzzy curves.
"""

import itertools
import functools
from collections import Counter, namedtuple
from collections.abc import Sized
import logging

from quicktions import Fraction

from .utils import get_lcm, get_int_cube_with_cache, get_int_time_with_cache, gen_faces
from . import _sat_adapters
from .curves import Curve, CurvePoint
from .subsets import Point
from ._bounded_items_heap import BoundedItemsHeap


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
        self._hash = hash(self._data())

    def _data(self):
        return self.cnums, self.cubes, self.div, self.dim  # dim is least used in comparison

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return self._hash

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

        orig_cnum = prev_spec.base_map.apply_cnum(self.curve.genus, self.pos.cnums[-1])  # cube in orig curve to divide

        for orig_spec in self.curve.gen_allowed_specs(active_pnum, orig_cnum):
            specified_curve = self.curve.specify(active_pnum, orig_cnum, orig_spec)
            # first we transform proto in orig, then apply prev map to the whole picture
            last_curve_proto = (prev_spec.base_map * orig_spec.base_map) * self.curve.patterns[orig_spec.pnum].proto
            for cnum, cube in enumerate(last_curve_proto):
                yield _CurvePiece(specified_curve, self.pnum, self.pos.specify(cnum, cube))


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
        # we do not look at junc depth as we are not interested in subdivisions of the curve only
        # instead we just analyze all its junctions
        return min(self.pos1.depth, self.pos2.depth)

    def divide_balanced(self):
        # Divide one of the fractions keeping pair balanced: depth1 == depth2 or depth1 == depth2 + 1.

        # use curve from divided piece because it has specified curve
        if self.pos1.depth > self.pos2.depth:
            for subpiece in self._get_piece(2).divide():
                yield _CurvePiecePair(subpiece.curve, self.junc, self.pos1, subpiece.pos)
        else:
            for subpiece in self._get_piece(1).divide():
                yield _CurvePiecePair(subpiece.curve, self.junc, subpiece.pos, self.pos2)


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

    def __init__(self, ratio_func, cache_max_size=2**18):
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
        self._get_pos_bounds = functools.lru_cache(cache_max_size)(self._get_pos_bounds)  # TODO: use methodtools?

    def _get_bounds(self, pair, curve_points=None):
        # Get lower and upper bounds for max ratio of given fractions pair:
        #   WD(f1,f2) := sup ||gamma(s)-gamma(t)||^d/|s-t|:  gamma(s) in f1, gamma(t) in f2
        # curve_points: dict {pnum: curve points list}
        # Returns triple (lo, up, argmax), argmax only if curve_points is set

        pos1, pos2 = pair.pos1, pair.pos2
        self.sum_stats['pair_depth'] += pair.depth
        self.max_stats['pair_depth'] = max(self.max_stats.get('pair_depth', 0), pair.depth)

        if curve_points is not None:
            # to avoid using curve in cached _get_pos_bounds method, we rotate points here
            sp1, sp2 = pair.get_last_specs()
            pts1 = tuple(sp1.base_map * pt for pt in curve_points[sp1.pnum])
            pts2 = tuple(sp2.base_map * pt for pt in curve_points[sp2.pnum])
        else:
            pts1 = pts2 = None

        return self._get_pos_bounds(pair.junc, pos1, pos2, pts1, pts2)

    def _get_pos_bounds(self, junc, pos1, pos2, pts1, pts2):
        dim, N = pos1.dim, pos1.div
        use_curve_points = (pts1 is not None)

        # these are integer positions in original curve patterns
        # we will transform them to absolute coords
        x1, t1 = pos1.get_int_coords()
        x2, t2 = pos2.get_int_coords()

        # scales:
        div1 = N**pos1.depth
        genus1 = div1**dim
        div2 = N**pos2.depth
        genus2 = div2**dim

        # junc: apply base_maps to coordinates
        jbm1, jbm2 = junc.spec1.base_map, junc.spec2.base_map
        x1, t1 = jbm1.apply_cube(div1, x1), jbm1.apply_cnum(genus1, t1)
        x2, t2 = jbm2.apply_cube(div2, x2), jbm2.apply_cnum(genus2, t2)
        if use_curve_points:
            pts1 = [jbm1 * pt for pt in pts1]
            pts2 = [jbm2 * pt for pt in pts2]

        # common scale
        if pos1.depth == pos2.depth:
            mx2 = mt2 = 1
        elif pos1.depth == pos2.depth + 1:
            mx2, mt2 = N, N**dim
            x2 = [xj * mx2 for xj in x2]
            t2 *= mt2
            if use_curve_points:
                pts2 = [pt.scale(mx2) for pt in pts2]
        else:
            raise ValueError("Unbalanced positions!")

        mx, mt = div1, genus1

        # now we have the following integer coordinates:
        # cube1: x1j <= xj <= x1j + 1    -- [0,  1]^d cube inside [0, mx]^d
        # cube2: x2j <= xj <= x2j + mx2  -- [0,mx2]^d cube inside [0, mx]^d, + shift junc_dx * mx
        #
        # time1: t1 <= t <= t1 + 1
        # time2: t2 <= t <= t2 + mt2,  + shift junc_dt * mt

        # junc: shifts
        x2 = [x2j + dxj * mx for x2j, dxj in zip(x2, junc.delta_x)]
        t2 += junc.delta_t * mt

        max_dx = [max(abs(x1j - x2j + 1), abs(x1j - x2j - mx2)) for x1j, x2j in zip(x1, x2)]

        max_dt = t2 + mt2 - t1  # max(t_2 - t_1)
        min_dt = t2 - (t1 + 1)  # min(t_2 - t_1)

        # also ratio(min_dx, min_dt) is a lower bound, but inefficient
        lo = self.ratio_func(dim, max_dx, max_dt)
        up = self.ratio_func(dim, max_dx, min_dt)

        argmax = None
        if use_curve_points:
            for (x1rel, t1rel), (x2rel, t2rel) in itertools.product(pts1, pts2):
                x1_point = [x1j + x1relj for x1j, x1relj in zip(x1, x1rel)]
                t1_point = t1 + t1rel
                x2_point = [x2j + x2relj for x2j, x2relj in zip(x2, x2rel)]
                t2_point = t2 + t2rel

                dx = [x1j - x2j for x1j, x2j in zip(x1_point, x2_point)]
                dt = t2_point - t1_point
                lo_point = self.ratio_func(dim, dx, dt)

                if lo_point > lo or argmax is None:
                    lo = lo_point
                    x1_real = [Fraction(x1j, mx) for x1j in x1_point]
                    t1_real = Fraction(t1_point, mt)
                    x2_real = [Fraction(x2j, mx) for x2j in x2_point]
                    t2_real = Fraction(t2_point, mt)
                    argmax = {'x1': x1_real, 't1': t1_real, 'x2': x2_real, 't2': t2_real, 'junc': junc, 'pos1': pos1, 'pos2': pos2}

        return lo, up, argmax

    def get_info(self):
        return {
            'bounds_cache': self._get_pos_bounds.cache_info(),
            'sum_stats': self.sum_stats,
            'max_stats': self.max_stats,
        }

    def _create_tree(self, curve, good_threshold=None, bad_threshold=None, **tree_kwargs):
        # Create initial pairs from a curve: for all junctions we take pairs of non-adjacent fractions
        G = curve.genus
        tree = BoundedItemsHeap(**tree_kwargs)
        tree.set_good_threshold(good_threshold)
        tree.set_bad_threshold(bad_threshold)

        pair_data = []
        for junc in curve.gen_auto_junctions():
            for cnum1 in range(G):
                for cnum2 in range(cnum1 + 2, G):
                    pair_data.append((junc, cnum1, cnum2))

        for junc in curve.gen_regular_junctions():
            last_cnum1 = 0 if junc.spec1.base_map.time_rev else G - 1
            first_cnum2 = G - 1 if junc.spec2.base_map.time_rev else 0
            for cnum1 in range(G):
                for cnum2 in range(G):
                    if (cnum1, cnum2) != (last_cnum1, first_cnum2):
                        pair_data.append((junc, cnum1, cnum2))

        for junc, cnum1, cnum2 in pair_data:
            pair = _CurvePiecePair.init_first_order(curve, junc, cnum1, cnum2)
            self._push_tree(tree, pair)

        return tree

    # items for heap = pairs of curve fractions ( _CurvePiecePair)
    # together with lo/up bounds on their dilation
    _BoundedPair = namedtuple('_BoundedPair', ['pair', 'lo', 'up', 'argmax'])

    def _push_tree(self, tree, pair):
        lo, up, argmax = self._get_bounds(pair, curve_points=tree.stash)
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

    def estimate_dilation_regular(self, curve, rel_tol_inv=100, max_iter=None, use_face_moments=False, face_dim=0, max_depth=None):
        """
        Estimate dilation for a regular peano curve (class Curve).

        Args:
            curve: Curve instance, fully defined polyfractal curve
            rel_tol_inv: inverted relative tolerance (may be set to None)
            max_iter: limit for subdivisions iters
            use_face_moments: use moments for dilation lower bounds (first+last)
            face_dim: dimension of faces for moments
            max_depth: allow not to consider pairs of higher depth
              note that fraction depth = min(position1 depth, position2 depth); junc depth is not counted!
              if dilation is attained at pair of fractions of depth <= max_depth,
              then resulting "lo" will be exact

        Returns:
            dict with keys:
            'lo': lower bound - dilation(curve) >= lo
            'up': upper bound - dilation(curve) <= up
            'argmax': pair of points where lo is achieved (if use_vertex_moments is set)
        """

        tree_kwargs = {'keep_max_lo_item': True}
        if use_face_moments:
            pts = {}
            for pnum in range(curve.mult):
                face_pts_set = set()
                for face in gen_faces(curve.dim, face_dim):
                    for last in [False, True]:
                        face_pts_set.add(curve.get_face_moment(face, pnum=pnum, last=last, find_point=True))
                pts[pnum] = tuple(sorted(face_pts_set))
            tree_kwargs['stash'] = pts

        if rel_tol_inv is not None:
            tolerance = Fraction(rel_tol_inv + 1, rel_tol_inv)

        # This is a basic algorithm that does not require SAT-solvers.
        # We maintain a "tree" (BoundedItemsHeap instance) of pairs (_BoundedPair)
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
                pairs_tree.cleanup()
                depth = min(item.pair.depth for item in pairs_tree.active_items())
                if depth > max_depth:
                    break

            self._divide_tree(pairs_tree)

            item = pairs_tree.max_lo_item
            if item.lo > curr_lo:
                curr_lo, argmax = item.lo, item.argmax
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
                              max_iter=None, sat_iter_multiplier=1.3,
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
            good_threshold: fraction, good curves are those with dilation <= good_thr
            bad_threshold: fraction, bad curves are those with dilation >= bad_thr
              It is required that bad_threshold < good_threshold
            max_iter: max number of divisions; raise exception if maximum iterations reached
            sat_iter_multiplier: X; we call sat solver on every X**k iteration
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
        # We grow the pairs tree with fixed good and bad thresholds,
        # all bad pairs are added to the list of forbidden configurations,
        # propagated into a list of boolean clauses in sat adapter.
        #
        # If we can't find a model, then all curves are necessarily bad
        # (i.e., all possible configurations are forbidden)
        # If there are no active pairs in the tree and adapter finds a model,
        # then the corresponding curve is good because it avoids bad pairs
        # and all other curve's pairs are good

        thrs = {'good_threshold': good_threshold, 'bad_threshold': bad_threshold}
        if init_pairs_tree is None:
            pairs_tree = self._create_tree(curve, **thrs)
        else:
            pairs_tree = init_pairs_tree.copy_and_cleanup(**thrs)

        if init_sat_adapter is None:
            adapter = _sat_adapters.CurveSATAdapter(curve=curve)
        else:
            adapter = _sat_adapters.CurveSATAdapter(adapter=init_sat_adapter)

        self._forbid(pairs_tree, adapter)

        no_model = None
        iter_no = 0
        sat_iter = 1
        while pairs_tree.has_items():
            iter_no += 1
            if (max_iter is not None) and iter_no > max_iter:
                raise self._RunOutOfIterationsException()

            self._divide_tree(pairs_tree)
            self._forbid(pairs_tree, adapter)
            self.sum_stats['divide_iter'] += 1

            if iter_no < sat_iter:
                continue

            # here we try sat solver
            sat_iter = int(sat_iter * sat_iter_multiplier) + 1
            logger.debug('iter %d, try SAT solver; stats: %s, %s', iter_no, self.sum_stats, self.max_stats)
            self.sum_stats['sat.solve_calls'] += 1
            if not adapter.solve():
                no_model = True
                break

        self._add_tree_stats(pairs_tree)
        self._update_max_stats('sat.problem_size', adapter.get_problem_size())

        if no_model or not adapter.solve():
            return None

        return adapter.get_model_curve()

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
              all bounds are Fraction instances
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
              increase to reduce mem usage
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
            epoch += 1
            curr_rel_tol_inv *= rel_tol_inv_mult
            total = len(active) if isinstance(active, list) else -1
            new_active = []  # heap of CurveItem
            for cnt, item in enumerate(active):
                if epoch == 1:
                    self.sum_stats['seen_pcurve'] += 1
                    self.sum_stats['possible_curves'] += item.curve.count_possible_curves()
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
