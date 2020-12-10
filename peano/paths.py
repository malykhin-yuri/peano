from collections import defaultdict
import logging
import itertools

from sympy import Rational

from .base_maps import BaseMap
from .subsets import Gate, Point, Portal
from .node_paths import NodePathTree
from .utils import combinations_product


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
        cubes = [base_map.apply_cube(self.div, cube) if cube is not None else None for cube in self]
        if base_map.time_rev:
            cubes = reversed(cubes)
        return type(self)(self.dim, self.div, cubes)

    @classmethod
    def parse_basis(cls, chain_code):
        """
        Convert chain code like 'ijK' to curve prototype.

        We assume that we start from zero cube!  TODO : do not assume:)
        """
        dim = len(set(''.join(chain_code).lower()))

        assert dim <= 6
        letters = 'ijklmn'

        vect_dict = {}
        for k in range(dim):
            coord = [0]*dim
            coord[k] = 1
            vect_dict[letters[k]] = coord
            vect_dict[letters[k].upper()] = [-m for m in coord]

        def diag_coord(vector):
            arg = [vect_dict[k] for k in vector]
            coord = list(map(sum,zip(*arg)))
            return coord

        proto = [list(map(vect_dict.get,chain_code)) if len(chain_code) == 1 else diag_coord(m) for m in chain_code]

        proto = [[0] * dim] + proto
        for l in range(len(proto)-1):
            proto[l+1] = [c + d for c, d in zip(proto[l], proto[l+1])]

        div = 1 + max(cj for cube in proto for cj in cube)

        return cls(dim, div, proto)


class PortalPath:
    """Prototype + portals."""

    def __init__(self, proto, portals):
        self.proto = proto
        self.dim = proto.dim
        self.div = proto.div
        self.portals = tuple(portals)

        entr = portals[0].entrance.map_to_cube(self.div, proto[0])
        exit = portals[-1].exit.map_to_cube(self.div, proto[-1])
        self.portal = Portal(entr, exit)

    def __rmul__(self, base_map):
        new_portals = [base_map * portal for portal in self.portals]
        if base_map.time_rev:
            new_portals.reverse()
        return type(self)(base_map * self.proto, new_portals)

    def reversed(self):
        return BaseMap.id_map(self.dim).reversed_time() * self

    def _data(self):
        return self.portal, self.proto, self.portals

    def __lt__(self, other):
        return self._data() < other._data()

    def __eq__(self, other):
        return self._data() == other._data()

    def __hash__(self):
        return hash(self._data())


class CurvePath(PortalPath):
    """Prototype + gates. Legacy. TODO: решить, оставляем или выпиливаем"""

    def __init__(self, proto, gates):
        super().__init__(proto, gates)
        self.gate = Gate(self.portal.entrance, self.portal.exit)
        self.gates = self.portals


class PathsGenerator:
    """Generate paths with given portals.

    Given a list of portals (pairs subsets of [0,1]^d),
    generate PortalPaths, i.e. paths such that portal exit
    in each fraction intersects (or equals) portal entrance in next fraction.
    """

    def __init__(self, dim, div, portals=None, hdist=None, max_cdist=None, mode='auto'):
        """
        Init paths generator.

        dim, div    --  subj
        portals     --  list of portals
        hist        --  one gate: (0,..,0) -> (0,..,0,1,1,..,1) with k ones
        max_cdist   --  maximum l1-distance between cubes
        mode        --  'auto'|'intersects'|'equals' - condition on portals
        TODO: mode Не используется, выпилим?
        """

        self.dim = dim
        self.div = div

        if portals is None:
            entrance = Point((Rational(0, 1),) * dim)
            exit = Point((Rational(0, 1),) * (dim - hdist) + (Rational(1, 1),) * hdist)
            portals = [Gate(entrance, exit)]
        self.portals = portals

        if mode == 'auto':
            if all(all(isinstance(subset, Point) for subset in [portal.entrance, portal.exit]) for portal in portals):
                mode = 'equals'
            else:
                mode = 'intersects'
        self.mode = mode

        entr2exits = defaultdict(dict)
        for bm in BaseMap.gen_base_maps(dim, time_rev=False):
            for portal in portals:
                bm_entr = bm * portal.entrance
                bm_exit = bm * portal.exit
                entr2exits[bm_entr][bm_exit] = 1
                entr2exits[bm_exit][bm_entr] = 1
        self.entr2exits = entr2exits
        self.next_dict = self.get_next_dict(max_cdist)

    def gen_intersected(self, subset):
        """
        Given a subset in [0,1]^d, find portals entrances that intersects it.
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
        Get a dict: exit_subset => [(cube_delta, new_portal), ...]

        Params:
            max_cdist:  do not allow cube changes greater than it
        """

        # there may be too many exit sets, work with standard
        std_exits = set()
        for portal in self.portals:
            std_exits.add(portal.exit.std())
            std_exits.add(portal.entrance.std())  # "no time_rev" is not supported

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
                        next_pos.append((next_cube, Portal(next_entr, next_exit)))

                result[exit] = next_pos

        return result

    def get_paths_example(self, parents=None, **kwargs):
        """Generate one paths tuple."""
        examples = []
        if parents is None:
            parents = (None,) * len(self.portals)
        for portal, parent in zip(self.portals, parents):
            path = next(self.generate_paths_generic(portal=portal, parent=parent, **kwargs), None)
            if path is None:
                return None
            examples.append(path)
        return examples

    def generate_paths(self, parents=None, std=False, **kwargs):
        """
        Generate tuples of paths (p_1,..,p_g) for self portals.

        parents  --  additional restriction for paths (see generic method)
        """
        portals = self.portals
        if parents is None:
            parents = (None,) * len(portals)
        kwargs['std'] = std

        # only for one portal we do not consume path into memory
        if len(portals) == 1:
            for path in self.generate_paths_generic(portals[0], parent=parents[0], **kwargs):
                yield (path,)
            return

        paths_dict = {}
        restrictions = list(zip(portals, parents))
        for portal, parent in set(restrictions):
            paths_dict[portal, parent] = list(self.generate_paths_generic(portal=portal, parent=parent, **kwargs))

        logging.info('generate_paths counts: %s', [len(paths_dict[r]) for r in restrictions])

        if std:
            yield from combinations_product(restrictions, paths_dict)
        else:
            path_lists = [paths_dict[r] for r in restrictions]
            yield from itertools.product(*path_lists)


    def generate_paths_generic(self, portal=None, parent=None, std=False, **kwargs):
        """
        Generate PortalPath with given (global) portal using all self.portals.

        This method is usually called for portal in self.portals, but this is not required.
        parent -- PortalPath such that generated paths must be consistent with it (same proto & intersecting portals)
        std -- try to standartize path (minimize keeping parent/portal)
        """
        # will use NodePathTree with state = portal
        start = []
        finish = []

        if parent is None:
            #assert parent is None  # TODO fix
            def gen_next(node):
                return self.next_dict[node.state.exit]

            def gen_prev(node):
                for delta, portal in self.next_dict[node.state.entrance]:
                    yield delta, ~portal

            tree = NodePathTree(dim=self.dim, div=self.div, next_func=gen_next, prev_func=gen_prev)

            for cube, cube_subset in portal.entrance.divide(self.div):
                for start_entr in self.gen_intersected(cube_subset):
                    for start_exit in self.entr2exits[start_entr]:
                        logging.debug('start at %s: %s -> %s', cube, start_entr, start_exit)
                        start.append(tree.init_path(cube, state=Portal(start_entr, start_exit)))

            for cube, cube_subset in portal.exit.divide(self.div):
                for finish_entr in self.gen_intersected(cube_subset):
                    for finish_exit in self.entr2exits[finish_entr]:
                        logging.debug('finish at %s: %s -> %s', cube, finish_exit, finish_entr)
                        finish.append(tree.init_path(cube, Portal(finish_exit, finish_entr)))
        else:
            def check_parent(cnum, cube, portal):
                """Check that portal at given cube is consistent with parent"""
                return cube == parent.proto[cnum] and portal.intersects(parent.portals[cnum])

            cube2cnum = {cube: cnum for cnum, cube in enumerate(parent.proto)}
            def gen_next(node):
                next_cnum = cube2cnum[node.cube] + 1  # currently we are at correct cube
                for delta, next_portal in self.next_dict[node.state.exit]:
                    next_cube = tuple(cj + dj for cj, dj in zip(node.cube, delta))
                    if check_parent(next_cnum, next_cube, next_portal):
                        yield delta, next_portal

            def gen_prev(node):
                next_cnum = cube2cnum[node.cube] - 1
                for delta, next_portal in self.next_dict[node.state.entrance]:
                    next_cube = tuple(cj + dj for cj, dj in zip(node.cube, delta))
                    if check_parent(next_cnum, next_cube, ~next_portal):
                        yield delta, ~next_portal

            tree = NodePathTree(dim=self.dim, div=self.div, next_func=gen_next, prev_func=gen_prev)

            start_cube = parent.proto[0]
            for start_entr in self.gen_intersected(parent.portals[0].entrance):
                for start_exit in self.entr2exits[start_entr]:
                    start_portal = Portal(start_entr, start_exit)
                    if check_parent(0, start_cube, start_portal):
                        start.append(tree.init_path(start_cube, start_portal))

            finish_cube = parent.proto[-1]
            for finish_entr in self.gen_intersected(parent.portals[-1].exit):
                for finish_exit in self.entr2exits[finish_entr]:
                    finish_portal = Portal(finish_exit, finish_entr)
                    if check_parent(-1, finish_cube, finish_portal):
                        finish.append(tree.init_path(finish_cube, finish_portal))

        if std:
            seen_paths = set()
            bms = BaseMap.gen_base_maps(self.dim)
            keep = parent if parent is not None else portal
            std_bms = [bm for bm in bms if bm * keep == keep]

        for path_data in tree.grow(start, ends=finish, **kwargs):
            result_path = PortalPath(
                proto=Proto(self.dim, self.div, [head.cube for head in path_data]),
                portals=[head.state for head in path_data],
            )
            if std:
                std_path = min(bm * result_path for bm in std_bms)
                if std_path not in seen_paths:
                    yield std_path
                    seen_paths.add(std_path)
                continue
            else:
                yield result_path
