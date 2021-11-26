from collections import defaultdict
import logging
import itertools
import re

from quicktions import Fraction

from .base_maps import BaseMap
from .subsets import Link, Point
from ._cube_path_trees import CubePathTree
from .utils import combinations_product, BASIS_LETTERS


logger = logging.getLogger(__name__)


class Proto(tuple):
    """
    Curve prototype -- sequence of cubes.

    We allow None-s for some of cubes, to support usage of get_entrance/get_exit methods.
    """

    def __new__(cls, dim, div, cubes):
        cubes = tuple(tuple(cube) if cube is not None else None for cube in cubes)
        obj = super().__new__(cls, cubes)
        obj.dim = dim
        obj.div = div
        return obj

    def __rmul__(self, base_map):
        """Apply a base map to prototype."""
        if base_map.dim != self.dim:
            raise Exception("Incompatible base map!")
        src_cubes = reversed(self) if base_map.time_rev else self
        cubes = (base_map.apply_cube(self.div, cube) if cube is not None else None for cube in src_cubes)
        return Proto(self.dim, self.div, cubes)

    def __invert__(self):
        """Time-reversed prototype."""
        return Proto(self.dim, self.div, reversed(self))

    @classmethod
    def parse(cls, chain_code):
        """
        Convert chain code like 'ijK' to curve prototype.

        If proto is not facet-continuous, use grouping: i(jk)J, here (jk) means j+k
        """
        chain_groups = [grp.strip('()') for grp in re.findall('\w|\(\w+\)', chain_code)]
        dim = len(set(c.lower() for grp in chain_groups for c in grp))
        assert dim <= len(BASIS_LETTERS)
        l2i = {l: i for i, l in enumerate(BASIS_LETTERS)}

        cubes = [(0,) * dim]
        for bases in chain_groups:
            delta = [0] * dim
            for basis in bases:
                i = l2i[basis.lower()]
                delta[i] = 1 if basis.islower() else -1

            next_cube = tuple(cj + dj for cj, dj in zip(cubes[-1], delta))
            cubes.append(next_cube)

        # fix if we start not from zero cube
        min_cube = min(cube for cube in cubes)
        cubes = [tuple(cj - mj for cj, mj in zip(cube, min_cube)) for cube in cubes]

        div = 1 + max(cj for cube in cubes for cj in cube)
        return cls(dim, div, cubes)

    def __str__(self):
        res = []
        for idx in range(len(self) - 1):
            delta = (nj - cj for nj, cj in zip(self[idx+1], self[idx]))
            bases = []
            for x, letter in zip(delta, BASIS_LETTERS):
                if x == 1:
                    bases.append(letter)
                elif x == -1:
                    bases.append(letter.upper())
            res.append(bases[0] if len(bases) == 1 else '({})'.format(''.join(bases)))
        return ''.join(res)


class Path:
    """
    Prototype with links.

    In each fraction, link defines entrance and exit subsets (relative to the fraction).
    Most important is the case of point links; such path may be also called
    "Pointed prototype"; it is the basis for dilation estimation using SAT-solvers.
    """
    def __init__(self, proto, links):
        self.proto = proto
        self.dim = proto.dim
        self.div = proto.div
        self.links = tuple(links)

        entr = self.links[0].entrance.map_to_cube(self.div, proto[0])
        exit = self.links[-1].exit.map_to_cube(self.div, proto[-1])
        self.link = Link(entr, exit)

    def __rmul__(self, base_map):
        src_links = reversed(self.links) if base_map.time_rev else self.links
        new_links = (base_map * link for link in src_links)
        return Path(base_map * self.proto, new_links)

    def __invert__(self):
        return ~BaseMap.id_map(self.dim) * self

    def _data(self):
        return self.link, self.proto, self.links

    def __lt__(self, other):
        return self._data() < other._data()

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())

    def is_continuous(self):
        """
        Check if path is continuous.

        For pointed paths the notion of continuity is obvious.
        For generic paths we mean that in each fraction, exit must
        intersect entrance of the next fraction. Only such paths
        may "contain" a continuous curve.
        """
        prev_cube = self.proto[0]
        prev_link = self.links[0]
        for cube, link in zip(self.proto[1:], self.links[1:]):
            # check prev_link.exit ~~ link.entrance
            shift = [cj - pj for cj, pj in zip(cube, prev_cube)]
            if not prev_link.exit.intersects(link.entrance.transform(shift=shift)):
                return False
            prev_cube, prev_link = cube, link
        return True


class PathsGenerator:
    """Generate continuous paths with given links."""

    def __init__(self, dim, div, links, max_cdist=None):
        """
        Init paths generator.

        Args:
            dim, div: subj
            links: list of Link instances, used in two ways:
              * as a set of possible links for each fraction in generate_paths_generic
                note that here we use non-oriented links, i.e. always add their reverse
              * as a default list of global links in generate_paths method
            max_cdist: maximum l1-distance between adjacent cubes (limitation on prototypes)
        """
        self.dim = dim
        self.div = div
        self.links = tuple(links)

        self._init_entr2links()
        self._is_pointed = all(isinstance(entr, Point) for entr in self._entr2links)
        self._init_exit2next(max_cdist)

    def _init_entr2links(self):
        entr2exits = defaultdict(dict)  # to keep order
        for bm in BaseMap.gen_base_maps(self.dim, time_rev=False):  # optimization: do not use time_rev
            for link in self.links:
                bm_entr = bm * link.entrance
                bm_exit = bm * link.exit
                entr2exits[bm_entr][bm_exit] = 1
                entr2exits[bm_exit][bm_entr] = 1

        self._entr2links = {}
        for entr, exit_dict in entr2exits.items():
            self._entr2links[entr] = tuple(Link(entr, exit) for exit in exit_dict)

    def _gen_intersected_entr(self, subset):
        # Given a subset in [0,1]^d, find entrances that intersect it
        if self._is_pointed and isinstance(subset, Point):
            if subset in self._entr2links:
                yield subset
        else:
            for entr in self._entr2links:
                if entr.intersects(subset):
                    yield entr

    def _gen_intersected_links(self, subset):
        # Given a subset in [0,1]^d, find links with entrance that intersects it
        for entr in self._gen_intersected_entr(subset):
            yield from self._entr2links[entr]

    def _init_exit2next(self, max_cdist=None):
        # setup self._exit2next: exit_subset => [(cube_delta, new_link), ...]

        # there may be too many exit sets, work with standard
        std_exits = set()
        for link in self.links:
            std_exits.add(link.exit.std())
            std_exits.add(link.entrance.std())

        self._exit2next = {}
        bms = list(BaseMap.gen_base_maps(self.dim, time_rev=False))
        for std_exit in std_exits:
            std_next = []
            for cube, cube_subset in std_exit.gen_neighbours():
                if (max_cdist is not None) and sum(abs(cj) for cj in cube) > max_cdist:
                    continue
                for entr in self._gen_intersected_entr(cube_subset):
                    std_next.append((cube, entr))

            for bm in bms:
                exit = bm * std_exit
                if exit in self._exit2next:
                    continue

                next_pos = []
                for cube, entr in std_next:
                    bm_cube = bm.apply_cube_start(cube, 1)  # this does not change cdist!
                    for next_link in self._entr2links[bm * entr]:  # do not use _gen_intersected_links, optimization!
                        next_pos.append((bm_cube, next_link))

                self._exit2next[exit] = tuple(next_pos)

    def _get_restrictions(self, links, parents):
        if links is None and parents is None:
            links = self.links
            parents = (None,) * len(links)
        elif links is None:
            links = (None,) * len(parents)
        elif parents is None:
            parents = (None,) * len(links)
        else:
            raise ValueError("provide either links or parents!")
        return tuple((link, path) for link, path in zip(links, parents))

    def get_paths_example(self, links=None, parents=None, **kwargs):
        """
        Generate one paths tuple.

        Args:
            see generate_paths method

        Returns:
            tuple of paths if found else None
        """
        examples = []
        for link, parent in self._get_restrictions(links, parents):
            path = next(self.generate_paths_generic(link=link, parent=parent, **kwargs), None)
            if path is None:
                return None
            examples.append(path)
        return tuple(examples)

    def generate_paths(self, links=None, parents=None, **kwargs):
        """
        Generate tuples of paths (p_1,..,p_g) with given restrictions, using self.links in fractions.

        Args:
            links: tuple of links for each pattern
            parents: tuple of parents for each pattern (see generate_paths_generic)
              use self.links if not links nor patterns given
            **kwargs: other kwargs passed to generate_paths_generic

        Yields:
            tuples of paths
        """
        restrictions = self._get_restrictions(links, parents)

        # only for one link we do not consume path into memory
        if len(restrictions) == 1:
            link, parent = restrictions[0]
            for path in self.generate_paths_generic(link=link, parent=parent, **kwargs):
                yield (path,)
            return

        paths_dict = {}
        for link, parent in set(restrictions):
            paths_dict[link, parent] = list(self.generate_paths_generic(link=link, parent=parent, **kwargs))

        logger.info('generate_paths counts: %s', [len(paths_dict[r]) for r in restrictions])

        if kwargs.get('std'):
            yield from combinations_product(restrictions, paths_dict)
        else:
            path_lists = [paths_dict[r] for r in restrictions]
            yield from itertools.product(*path_lists)

    def generate_paths_generic(self, link=None, parent=None, std=False, **kwargs):
        """
        Generate Path with given global restriction using all self.links in fractions.

        Args:
            link: global link for path
              usually one of self.links, but this is not required
            parent: path such that generated paths must be consistent with it (same proto & intersecting links)
              exactly one of (link, parent) must be defined
            std: standartize paths (minimize using admissible base maps)
              with point self.links and point input link std guarantees to yield unique paths
            **kwargs: passed to tree.grow: start_max_count, finish_max_count

        Yields:
            paths (Path instances) that are "continuous", see Path.is_continuous method
        """
        if not ((link is None) ^ (parent is None)):
            raise ValueError("exactly one of (link, parent) must be defined!")

        # will use cube path tree with state = link

        if parent is not None:
            cube2cnum = {cube: cnum for cnum, cube in enumerate(parent.proto)}
            def check_parent(cnum, cube, link):
                return cube == parent.proto[cnum] and link.intersects(parent.links[cnum])

            def check_next(cube, delta, next_link, reverse=False):
                next_cnum = cube2cnum[cube] + (-1 if reverse else 1)
                next_cube = tuple(cj + dj for cj, dj in zip(cube, delta))
                return check_parent(next_cnum, next_cube, next_link)

        def gen_next(cube, link):
            for delta, next_link in self._exit2next[link.exit]:
                if (parent is not None) and (not check_next(cube, delta, next_link)):
                    continue
                yield delta, next_link

        def gen_prev(cube, link):
            for delta, next_link in self._exit2next[link.entrance]:
                prev_link = ~next_link
                if (parent is not None) and (not check_next(cube, delta, prev_link, reverse=True)):
                    continue
                yield delta, prev_link

        tree = CubePathTree(dim=self.dim, div=self.div, next_func=gen_next, prev_func=gen_prev)

        if parent is None:
            start_cube_sets = link.entrance.divide(self.div)
            finish_cube_sets = link.exit.divide(self.div)
        else:
            start_cube_sets = [(parent.proto[0], parent.links[0].entrance)]
            finish_cube_sets = [(parent.proto[-1], parent.links[-1].exit)]

        start = []
        for start_cube, start_set in start_cube_sets:
            for start_link in self._gen_intersected_links(start_set):
                if (parent is not None) and (not check_parent(0, start_cube, start_link)):
                    continue
                logger.debug('start at %s, %s', start_cube, start_link)
                start.append((start_cube, start_link))

        finish = []
        for finish_cube, finish_set in finish_cube_sets:
            for finish_link in self._gen_intersected_links(finish_set):
                finish_link = ~finish_link
                if (parent is not None) and (not check_parent(-1, finish_cube, finish_link)):
                    continue
                logger.debug('finish at %s, %s', finish_cube, finish_link)
                finish.append((finish_cube, finish_link))

        if std:
            seen_paths = set()
            bms = BaseMap.gen_base_maps(self.dim)
            keep = parent if parent is not None else link
            std_bms = [bm for bm in bms if bm * keep == keep]

        for path_data in tree.grow(start, finish, **kwargs):
            result_path = Path(
                proto=Proto(self.dim, self.div, [cube for cube, _ in path_data]),
                links=[state for _, state in path_data],
            )
            if std:
                std_path = min(bm * result_path for bm in std_bms)
                if std_path in seen_paths:
                    continue
                seen_paths.add(std_path)

            yield result_path
