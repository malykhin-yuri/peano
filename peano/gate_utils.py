from collections import Counter
import logging
import itertools

from .curves import FuzzyCurve
from .subsets import Link, HyperFaceDivSubset
from .paths import Proto, PathsGenerator
from .base_maps import BaseMap, Spec


class GatesGenerator:
    def __init__(self, dim, div, pcount, hyper=False):
        """
        Generate gates for given configuration, optimized for hypercurves.

        Args:
          dim, div, pcount  --  base configuration
          hyper  --  work with hypercurves only
        """
        self.dim = dim
        self.div = div
        self.pcount = pcount
        self.hyper = hyper
        self.stats = Counter()
        self.seen_gates = set()
        self.seen_std_gates = set()
        self.pnum_cnum_list = [(pnum, cnum) for pnum in range(pcount) for cnum in [0, -1]]

    def gen_gates(self, **kwargs):
        """
        Generate all gates such that there is at least one path with thems.

        For hyper, we consider only hypercurves (with entrance/exit strictly on hyperfaces).
        This class is optimized for hypercurves search.

        Otherwise, we consider non-internal curves, i.e., curves with all gates on cube boundary.
        Indeed, if curve in internal, then only internal curves can "use" it,
        and non-internal curves use each other and have lower ratio.

        All objects are standartized: obj -> std(obj) = min(bm * obj for bm in base_maps).
        """
        self.seen_gates.clear()
        self.seen_std_gates.clear()
        if self.hyper:
            yield from self._gen_hyper_gates(**kwargs)
        else:
            yield from self._gen_possible_gates(**kwargs)

    def _gen_hyper_gates(self, narrow_steps=1):
        """
        Generate gates for given configuration, optimized for hypercurves.

        Args:
          narrow_steps  --  TODO
        """
        dim, div, pcount = self.dim, self.div, self.pcount

        face0 = HyperFaceDivSubset(dim=dim, div=div, face=(0, 0))  # x0=0
        face1 = HyperFaceDivSubset(dim=dim, div=div, face=(0, 1))  # x0=1
        face2 = HyperFaceDivSubset(dim=dim, div=div, face=(1, 0))  # x1=0

        # we either go to the opposite face, either to the neighbour face, or return to the same face
        link_variants = [Link(face0, face0), Link(face0, face1), Link(face0, face2)]
        links_list = list(itertools.combinations_with_replacement(link_variants, r=pcount))

        for links_idx, links in enumerate(links_list):
            logging.info('processing global links %d of %d', links_idx + 1, len(links_list))
            pg = PathsGenerator(dim=dim, div=div, links=links)
            plist = list(pg.generate_paths(std=True))

            for narrow_idx in range(narrow_steps):
                new_plist = []
                for paths_idx, paths in enumerate(plist):
                    logging.info('processing narrow: %d of %d', paths_idx + 1, len(plist))
                    new_plist += self.get_narrow_paths(paths)
                logging.info('narrow step %d or %d: %d => %d', narrow_idx + 1, narrow_steps, len(plist), len(new_plist))
                plist = new_plist

            for paths_idx, paths in enumerate(plist):
                logging.info('processing paths %d of %d, stats: %s', paths_idx + 1, len(plist), self.stats)
                yield from self._gen_path_gates(paths)

    def _gen_possible_gates(self):
        dim, div, pcount = self.dim, self.div, self.pcount

        # curves are non-internal, so first and last cubes are on the boundary
        all_cubes = itertools.product(range(div), repeat=dim)
        face_cubes = (cube for cube in all_cubes if any(cj == 0 or cj == div - 1 for cj in cube))
        face_pairs = itertools.combinations(face_cubes, 2)

        def std_pair(pair):
            def gen_all():
                for bm in BaseMap.gen_base_maps(dim, time_rev=False):
                    c1 = bm.apply_cube(div, pair[0])
                    c2 = bm.apply_cube(div, pair[1])
                    yield (c1, c2)
                    yield (c2, c1)
            return min(gen_all())

        std_pairs = sorted(set(std_pair(pair) for pair in face_pairs))

        # as base_maps are not defined yet, we can permute pairs (i.e., assume that pairs tuple is sorted)
        # so we use combinations with replacement instead of product
        std_pairs_list = list(itertools.combinations_with_replacement(std_pairs, r=pcount))

        # to find gates, we should specify:
        # * first and last cubes in proto (for each pattern) -- will use std_pairs_list for that
        # * specs on that cubes -- see below

        def touch_same_hyperface(cube1, cube2):
            return any((c1j == c2j == 0) or (c1j == c2j == div-1) for c1j, c2j in zip(cube1, cube2))

        for pairs_idx, pairs in enumerate(std_pairs_list):
            logging.info('processing cube pairs: %d of %d', pairs_idx + 1, len(std_pairs_list))
            protos = []
            for pair in pairs:
                proto = [None] * (div**dim)
                proto[0], proto[-1] = pair
                protos.append(Proto(dim, div, proto))

            spec_dict = {}
            for pnum, cnum in self.pnum_cnum_list:
                specs = []
                for pn in range(pcount):
                    for bm in BaseMap.gen_base_maps(dim):
                        if touch_same_hyperface(protos[pnum][cnum], (bm * protos[pn])[cnum]):
                            spec = Spec(bm, pnum=pn)
                            specs.append(spec)
                spec_dict[pnum, cnum] = specs
            yield from self._check_variants(protos, spec_dict)

    def _gen_path_gates(self, paths):
        # to determine: specs for first and last for all patterns
        # тонкий момент: в списке path.links мы не знаем, какой link из какого шаблона
        # а вот path.link привязан именно к этому шаблону
        spec_dict = {}
        for pnum, cnum in self.pnum_cnum_list:
            specs = []
            for pn in range(self.pcount):
                for bm in paths[pn].link.argmul_intersect(paths[pnum].links[cnum]):
                    spec = Spec(bm, pnum=pn)
                    specs.append(spec)
            spec_dict[pnum, cnum] = specs

        yield from self._check_variants([path.proto for path in paths], spec_dict)

    def _check_variants(self, protos, spec_dict):
        variants = [spec_dict[pnum, cnum] for pnum, cnum in self.pnum_cnum_list]
        total = 1
        for v in variants:
            total *= len(v)
        for specs_idx, specs in enumerate(itertools.product(*variants)):
            if (specs_idx + 1) % 1000 == 0:
                logging.info('processing variants: %d of %d', specs_idx + 1, total)
            gate = self.check_and_std(protos, specs)
            if gate is not None:
                yield gate

    def check_and_std(self, protos, spec_list):
        dim, div, pcount = self.dim, self.div, self.pcount

        pattern_specs = [[None] * (div**dim) for _ in range(pcount)]
        for (pnum, cnum), spec in zip(self.pnum_cnum_list, spec_list):
            pattern_specs[pnum][cnum] = spec

        patterns = [(proto, specs) for proto, specs in zip(protos, pattern_specs)]
        curve = FuzzyCurve(dim=dim, div=div, patterns=patterns)

        gates = []
        for pnum in range(pcount):
            entr = curve.get_entrance(pnum)
            if (self.hyper and entr.face_dim() != dim-1) or (entr.face_dim() == dim):
                return

            exit = curve.get_exit(pnum)
            if (self.hyper and exit.face_dim() != dim-1) or (exit.face_dim() == dim):
                return

            gates.append(Link(entr, exit))

        gates = tuple(gates)
        if gates in self.seen_gates:
            return

        self.seen_gates.add(gates)
        self.stats['new_path_gate'] += 1

        std_gates = tuple(sorted(gate.std() for gate in gates))

        if std_gates in self.seen_std_gates:
            return
        self.seen_std_gates.add(std_gates)
        self.stats['new_std_gate'] += 1

        # check that there is at least one path with given gates
        pg = PathsGenerator(dim=dim, div=div, links=std_gates)

        # параметры -- дискуссионный вопрос
        if pg.get_paths_example(start_max_count=1000, finish_max_count=100000):
            logging.debug('GOOD gates: %s', [str(g) for g in std_gates])
            self.stats['new_good_gate'] += 1
            return std_gates
        else:
            logging.debug('BAD gates: %s', [str(g) for g in std_gates])

    def get_narrow_paths(self, paths):
        div = paths[0].proto.div
        dim = paths[0].proto.dim

        narrow_links = []
        for path in paths:
            narrow_links.append(path.link)

        npg = PathsGenerator(dim=dim, div=div, links=narrow_links)
        if not npg.get_paths_example(parents=paths):
            return

        narrows = []
        for path, narrow_link in zip(paths, narrow_links):
            seen = set()
            res = []
            for narrow_path in npg.generate_paths_generic(link=path.link, parent=path):
                key = (narrow_path.link, narrow_path.links[0], narrow_path.links[-1])
                if key not in seen:
                    seen.add(key)
                    res.append(narrow_path)
            narrows.append(res)

        # TODO: а здесь conbinations product не ?? или просто генерить со списком парентов?
        yield from itertools.product(*narrows)
