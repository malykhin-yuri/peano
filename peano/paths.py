from collections import defaultdict
import logging
import itertools
import re

from quicktions import Fraction

from .base_maps import BaseMap
from .subsets import Point, Link
from ._cube_path_trees import CubePathTree
from .utils import combinations_product


logger = logging.getLogger(__name__)


class Proto(tuple):
    """
    Curve prototype -- sequence of cubes.

    We allow None-s for some of cubes, to support usage of get_entrance/get_exit methods.
    """

    basis_letters = 'ijklmn'

    def __new__(cls, dim, div, cubes):
        cubes = tuple(tuple(cube) if cube is not None else None for cube in cubes)
        obj = super().__new__(cls, cubes)
        obj.dim = dim
        obj.div = div
        return obj

    def __rmul__(self, base_map):
        cubes = [base_map.apply_cube(self.div, cube) if cube is not None else None for cube in self]
        if base_map.time_rev:
            cubes = reversed(cubes)
        return type(self)(self.dim, self.div, cubes)

    def __invert__(self):
        """Time-reversed prototype."""
        return type(self)(self.dim, self.div, tuple(reversed(self)))

    @classmethod
    def parse(cls, chain_code):
        """
        Convert chain code like 'ijK' to curve prototype.
        If proto is not facet-continuous, use groups: i(jk)J
        """
        chain_groups = [grp.strip('()') for grp in re.findall('\w|\(\w+\)', chain_code)]
        dim = len(set(c.lower() for grp in chain_groups for c in grp))
        assert dim <= len(cls.basis_letters)
        l2i = {l: i for i, l in enumerate(cls.basis_letters)}

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
            delta = tuple(nj - cj for nj, cj in zip(self[idx+1], self[idx]))
            bases = []
            for x, letter in zip(delta, self.basis_letters):
                if x == 1:
                    bases.append(letter)
                elif x == -1:
                    bases.append(letter.upper())
            if len(bases) == 1:
                res.append(bases[0])
            else:
                res.append('(' + ''.join(bases) + ')')
        return ''.join(res)


class Path:
    """Prototype + links."""
    def __init__(self, proto, links):
        self.proto = proto
        self.dim = proto.dim
        self.div = proto.div
        self.links = tuple(links)

        entr = links[0].entrance.map_to_cube(self.div, proto[0])
        exit = links[-1].exit.map_to_cube(self.div, proto[-1])
        self.link = Link(entr, exit)

    def __rmul__(self, base_map):
        new_links = [base_map * link for link in self.links]
        if base_map.time_rev:
            new_links.reverse()
        return type(self)(base_map * self.proto, new_links)

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

    def is_pointed(self):
        return all(link.is_pointed() for link in self.links)

    def is_continuous(self):
        if not self.is_pointed():
            raise TypeError("Continuity is defined only for pointed paths")
        prev_cube = self.proto[0]
        prev_link = self.links[0]
        for cube, link in zip(self.proto[1:], self.links[1:]):
            # check prev_link.exit ~~ link.entrance
            shift = [cj - pj for cj, pj in zip(cube, prev_cube)]
            if prev_link.exit != link.entrance.transform(shift=shift):
                return False
            prev_cube, prev_link = cube, link
        return True


class PathsGenerator:
    """Generate paths with given links.

    Given a list of links (pairs subsets of [0,1]^d),
    generate Paths, i.e. paths such that link exit
    in each fraction intersects (or equals) link entrance in next fraction.
    """

    def __init__(self, dim, div, links=None, hdist=None, max_cdist=None, mode='auto'):
        """
        Init paths generator.

        dim, div    --  subj
        links     --  list of links
        hist        --  one gate: (0,..,0) -> (0,..,0,1,1,..,1) with k ones
        max_cdist   --  maximum l1-distance between cubes
        mode        --  'auto'|'intersects'|'equals' - condition on links
        TODO: mode Не используется, выпилим?
        """

        self.dim = dim
        self.div = div

        if links is None:
            entrance = Point((Fraction(0),) * dim)
            exit = Point((Fraction(0),) * (dim - hdist) + (Fraction(1, 1),) * hdist)
            links = [Link(entrance, exit)]
        self.links = links

        if mode == 'auto':
            if all(link.is_pointed() for link in links):
                mode = 'equals'
            else:
                mode = 'intersects'
        self.mode = mode

        entr2exits = defaultdict(dict)
        for bm in BaseMap.gen_base_maps(dim, time_rev=False):
            for link in links:
                bm_entr = bm * link.entrance
                bm_exit = bm * link.exit
                entr2exits[bm_entr][bm_exit] = 1
                entr2exits[bm_exit][bm_entr] = 1
        self.entr2exits = entr2exits
        self.next_dict = self.get_next_dict(max_cdist)

    def gen_intersected(self, subset):
        """
        Given a subset in [0,1]^d, find links entrances that intersects it.
        """
        if self.mode == 'equals':
            if subset in self.entr2exits:
                yield subset
        else:
            for entr in self.entr2exits:
                if entr.intersects(subset):
                    yield entr

    def get_next_dict(self, max_cdist=None):
        """
        Get a dict: exit_subset => [(cube_delta, new_link), ...]

        Params:
            max_cdist:  do not allow cube changes greater than it
        """

        # there may be too many exit sets, work with standard
        std_exits = set()
        for link in self.links:
            std_exits.add(link.exit.std())
            std_exits.add(link.entrance.std())  # "no time_rev" is not supported

        result = {}
        bms = list(BaseMap.gen_base_maps(self.dim, time_rev=False))
        for std_exit in std_exits:
            std_next = []
            for cube, cube_subset in std_exit.gen_neighbours():
                if max_cdist is not None:
                    if sum(abs(cj) for cj in cube) > max_cdist:
                        continue
                for entr in self.gen_intersected(cube_subset):
                    std_next.append((cube, entr))

            for bm in bms:
                exit = bm * std_exit
                if exit in result:
                    continue

                next_pos = []
                for cube, entr in std_next:
                    next_cube = bm.apply_cube_start(cube, 1)  # this does not change cdist!
                    next_entr = bm * entr
                    for next_exit in self.entr2exits[next_entr]:
                        next_pos.append((next_cube, Link(next_entr, next_exit)))

                result[exit] = next_pos

        return result

    def get_paths_example(self, parents=None, **kwargs):
        """Generate one paths tuple."""
        examples = []
        if parents is None:
            parents = (None,) * len(self.links)
        for link, parent in zip(self.links, parents):
            path = next(self.generate_paths_generic(link=link, parent=parent, **kwargs), None)
            if path is None:
                return None
            examples.append(path)
        return examples

    def generate_paths(self, parents=None, std=False, **kwargs):
        """
        Generate tuples of paths (p_1,..,p_g) for self links.

        parents  --  additional restriction for paths (see generic method)
        """
        links = self.links
        if parents is None:
            parents = (None,) * len(links)
        kwargs['std'] = std

        # only for one link we do not consume path into memory
        if len(links) == 1:
            for path in self.generate_paths_generic(links[0], parent=parents[0], **kwargs):
                yield (path,)
            return

        paths_dict = {}
        restrictions = list(zip(links, parents))
        for link, parent in set(restrictions):
            paths_dict[link, parent] = list(self.generate_paths_generic(link=link, parent=parent, **kwargs))

        logger.info('generate_paths counts: %s', [len(paths_dict[r]) for r in restrictions])

        if std:
            yield from combinations_product(restrictions, paths_dict)
        else:
            path_lists = [paths_dict[r] for r in restrictions]
            yield from itertools.product(*path_lists)

    def generate_paths_generic(self, link=None, parent=None, std=False, **kwargs):
        """
        Generate Path with given (global) link using all self.links.

        This method is usually called for link in self.links, but this is not required.
        parent -- Path such that generated paths must be consistent with it (same proto & intersecting links)
        std -- try to standartize path (minimize keeping parent/link)
        """
        # will use CubePathTree with state = link
        start = []
        finish = []

        if parent is None:
            #assert parent is None  # TODO fix
            def gen_next(cube, link):
                return self.next_dict[link.exit]

            def gen_prev(cube, link):
                for delta, link in self.next_dict[link.entrance]:
                    yield delta, ~link

            tree = CubePathTree(dim=self.dim, div=self.div, next_func=gen_next, prev_func=gen_prev)

            for cube, cube_subset in link.entrance.divide(self.div):
                for start_entr in self.gen_intersected(cube_subset):
                    for start_exit in self.entr2exits[start_entr]:
                        logger.debug('start at %s: %s -> %s', cube, start_entr, start_exit)
                        start.append((cube, Link(start_entr, start_exit)))

            for cube, cube_subset in link.exit.divide(self.div):
                for finish_entr in self.gen_intersected(cube_subset):
                    for finish_exit in self.entr2exits[finish_entr]:
                        logger.debug('finish at %s: %s -> %s', cube, finish_exit, finish_entr)
                        finish.append((cube, Link(finish_exit, finish_entr)))
        else:
            def check_parent(cnum, cube, link):
                """Check that link at given cube is consistent with parent"""
                return cube == parent.proto[cnum] and link.intersects(parent.links[cnum])

            cube2cnum = {cube: cnum for cnum, cube in enumerate(parent.proto)}
            def gen_next(cube, link):
                next_cnum = cube2cnum[cube] + 1  # currently we are at correct cube
                for delta, next_link in self.next_dict[link.exit]:
                    next_cube = tuple(cj + dj for cj, dj in zip(cube, delta))
                    if check_parent(next_cnum, next_cube, next_link):
                        yield delta, next_link

            def gen_prev(cube, link):
                next_cnum = cube2cnum[cube] - 1
                for delta, next_link in self.next_dict[link.entrance]:
                    next_cube = tuple(cj + dj for cj, dj in zip(cube, delta))
                    if check_parent(next_cnum, next_cube, ~next_link):
                        yield delta, ~next_link

            tree = CubePathTree(dim=self.dim, div=self.div, next_func=gen_next, prev_func=gen_prev)

            start_cube = parent.proto[0]
            for start_entr in self.gen_intersected(parent.links[0].entrance):
                for start_exit in self.entr2exits[start_entr]:
                    start_link = Link(start_entr, start_exit)
                    if check_parent(0, start_cube, start_link):
                        start.append((start_cube, start_link))

            finish_cube = parent.proto[-1]
            for finish_entr in self.gen_intersected(parent.links[-1].exit):
                for finish_exit in self.entr2exits[finish_entr]:
                    finish_link = Link(finish_exit, finish_entr)
                    if check_parent(-1, finish_cube, finish_link):
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
                if std_path not in seen_paths:
                    yield std_path
                    seen_paths.add(std_path)
                continue
            else:
                yield result_path
