from collections import Counter
import logging
import itertools

from .curves import FuzzyCurve
from .subsets import Link, FacetDivSubset
from .paths import Proto, PathsGenerator
from .base_maps import BaseMap, Spec


logger = logging.getLogger(__name__)


class GatesGenerator:
    def __init__(self, dim, div, pcount, only_facet=False):
        """
        Generate gates for given configuration.

        Args:
            dim, div, pcount: base configuration
            only_facet: generate facet-gated curves only, i.e. having entrance/exit
              strictly on facets (hyperfaces)
        """
        self.dim = dim
        self.div = div
        self.pcount = pcount
        self.only_facet = only_facet
        self.stats = Counter()
        self._boundary_pnum_cnums = tuple((pnum, cnum) for pnum in range(pcount) for cnum in [0, -1])

    def gen_gates(self):
        """
        Generate all gates such that there is at least one path with them.

        Only non-internal curves are considered, i.e., curves with all gates on cube boundary.
        Indeed, if curve is internal, then only internal curves can "use" it,
        and non-internal curves use each other and have lower ratio.

        All objects are standartized: obj -> std(obj) = min(bm * obj for bm in base_maps).

        Yields:
            tuples of pointed Links instances
        """
        self._seen_gates = set()
        self._seen_std_gates = set()
        if self.only_facet:
            yield from self._gen_facet_gates()
        else:
            yield from self._gen_all_gates()

    def _gen_facet_gates(self):
        dim, div, pcount = self.dim, self.div, self.pcount

        facet0 = FacetDivSubset(dim=dim, div=div, facet=(0, 0))  # x0=0
        facet1 = FacetDivSubset(dim=dim, div=div, facet=(0, 1))  # x0=1
        facet2 = FacetDivSubset(dim=dim, div=div, facet=(1, 0))  # x1=0

        # we either go to the opposite facet, either to the neighbour facet, or return to the same facet
        link_variants = [Link(facet0, facet0), Link(facet0, facet1), Link(facet0, facet2)]
        links_list = list(itertools.combinations_with_replacement(link_variants, r=pcount))

        for links_idx, links in enumerate(links_list):
            logger.info('processing global links %d of %d', links_idx + 1, len(links_list))
            pg = PathsGenerator(dim=dim, div=div, links=links)
            plist = list(pg.generate_paths(std=True))

            # _get_narrow_paths uniqs by intermediate links, so we can't make 2+ narrowing steps!
            # this is optimized for 3d bifractals search; should experiment more ...
            for paths_idx, paths in enumerate(plist):
                logger.debug('processing paths %d of %d, stats: %s', paths_idx + 1, len(plist), self.stats)
                for narrow_paths in self._gen_narrow_paths(paths):
                    # to determine: specs for first and last for all patterns
                    specs_dict = {}
                    for pnum, cnum in self._boundary_pnum_cnums:
                        specs = []
                        for pn in range(self.pcount):
                            bms = narrow_paths[pn].link.argmul_intersect(narrow_paths[pnum].links[cnum])
                            specs += [Spec(bm, pnum=pn) for bm in bms]
                        specs_dict[pnum, cnum] = specs
                    yield from self._check_variants([path.proto for path in narrow_paths], specs_dict)

    def _gen_all_gates(self):
        dim, div, pcount = self.dim, self.div, self.pcount

        # curves are non-internal, so first and last cubes are on the boundary
        all_bms = list(BaseMap.gen_base_maps(dim=dim))
        all_cubes = itertools.product(range(div), repeat=dim)
        facet_cubes = (cube for cube in all_cubes if any(cj == 0 or cj == div - 1 for cj in cube))
        facet_pairs = itertools.combinations(facet_cubes, 2)

        def std_pair(pair):
            min_pair = pair
            for bm in BaseMap.gen_base_maps(dim, time_rev=False):
                c1 = bm.apply_cube(div, pair[0])
                c2 = bm.apply_cube(div, pair[1])
                min_pair = min(min_pair, (c1, c2), (c2, c1))
            return min_pair

        # as base_maps are not defined yet, we can standartize everything
        std_pairs = set(std_pair(pair) for pair in facet_pairs)
        std_pairs_list = list(itertools.combinations_with_replacement(sorted(std_pairs), r=pcount))

        # to find gates, we should specify:
        # * first and last cubes in proto (for each pattern) -- will use std_pairs_list for that
        # * specs on that cubes -- see below

        def touch_same_facet(cube1, cube2):
            return any((c1j == c2j == 0) or (c1j == c2j == div-1) for c1j, c2j in zip(cube1, cube2))

        for pairs_idx, pairs in enumerate(std_pairs_list):
            logger.info('processing cube pairs: %d of %d', pairs_idx + 1, len(std_pairs_list))
            protos = []
            for pair in pairs:
                proto = [None] * (div**dim)
                proto[0], proto[-1] = pair
                protos.append(Proto(dim, div, proto))

            specs_dict = {}
            for pnum, cnum in self._boundary_pnum_cnums:
                specs = []
                for pn in range(pcount):
                    specs += [Spec(bm, pnum=pn) for bm in all_bms if touch_same_facet(protos[pnum][cnum], (bm * protos[pn])[cnum])]
                specs_dict[pnum, cnum] = specs
            yield from self._check_variants(protos, specs_dict)

    def _check_variants(self, protos, specs_dict):
        variants = [specs_dict[pnum, cnum] for pnum, cnum in self._boundary_pnum_cnums]
        total = 1
        for v in variants:
            total *= len(v)
        for specs_idx, specs in enumerate(itertools.product(*variants)):
            if (specs_idx + 1) % 1000 == 0:
                logger.debug('processing variants: %d of %d', specs_idx + 1, total)
            gates = self._check_and_std(protos, specs)
            if gates is not None:
                yield gates

    def _check_and_std(self, protos, spec_list):
        dim, div, pcount = self.dim, self.div, self.pcount

        pattern_specs = [[None] * (div**dim) for _ in range(pcount)]
        for (pnum, cnum), spec in zip(self._boundary_pnum_cnums, spec_list):
            pattern_specs[pnum][cnum] = spec

        patterns = [(proto, specs) for proto, specs in zip(protos, pattern_specs)]
        curve = FuzzyCurve(dim=dim, div=div, patterns=patterns)

        gates = []
        for pnum in range(pcount):
            entr = curve.get_entrance(pnum)
            if (self.only_facet and entr.face_dim() != dim-1) or (entr.face_dim() == dim):
                return

            exit = curve.get_exit(pnum)
            if (self.only_facet and exit.face_dim() != dim-1) or (exit.face_dim() == dim):
                return

            gates.append(Link(entr, exit))

        gates = tuple(gates)
        if gates in self._seen_gates:
            return

        self._seen_gates.add(gates)
        self.stats['new_path_gate'] += 1

        std_gates = tuple(sorted(gate.std() for gate in gates))

        if std_gates in self._seen_std_gates:
            return
        self._seen_std_gates.add(std_gates)
        self.stats['new_std_gate'] += 1

        # check that there is at least one path with given gates
        # we can't use parent proto here because of gates "cache"
        pg = PathsGenerator(dim=dim, div=div, links=std_gates)

        # TODO: optimize params ?
        if not pg.get_paths_example(start_max_count=1000, finish_max_count=1000000):
            logger.debug('BAD gates: %s', [str(g) for g in std_gates])
            return

        logger.debug('GOOD gates: %s', [str(g) for g in std_gates])
        self.stats['new_good_gate'] += 1
        return std_gates

    def _gen_narrow_paths(self, paths):
        # paths = one paths tuple
        # narrowing idea: suppose we have a path with facet links;
        # this path's global link is more "narrow" - a div-subset of facet
        # so we output paths with that narrow links
        #
        # important optimization: we uniq paths by intermediate links, they are not used later

        dim, div = paths[0].proto.dim, paths[0].proto.div

        narrow_links = tuple(path.link for path in paths)
        npg = PathsGenerator(dim=dim, div=div, links=narrow_links)
        if not npg.get_paths_example(parents=paths):
            return

        narrows = []
        for path, narrow_link in zip(paths, narrow_links):
            seen = set()
            res = []
            for narrow_path in npg.generate_paths_generic(parent=path):
                key = (narrow_path.link, narrow_path.links[0], narrow_path.links[-1])
                if key not in seen:  # main optimization
                    seen.add(key)
                    res.append(narrow_path)
            narrows.append(res)

        yield from itertools.product(*narrows)
