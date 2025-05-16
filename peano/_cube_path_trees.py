from __future__ import annotations
import itertools
import logging
from collections import namedtuple, defaultdict
from collections.abc import Hashable, Set
from dataclasses import dataclass
from typing import Literal, Callable, Iterable


logger = logging.getLogger(__name__)

type CubeType = tuple[int, ...]

@dataclass(frozen=True)
class _CubeNode[S: Hashable]:
    # Node of the partial path in the cube lattice.
    cube: CubeType  # cube position
    state: S  # any "state", e.g. exit point in the cube

    def __str__(self) -> str:
        return '{} @{}'.format(self.cube, self.state)


@dataclass
class _CubePath[S]:
    # Path of nodes, with usual head+link_to_tail representation.
    head: _CubeNode[S]
    length: int
    prev: _CubePath[S] | None  # link to tail path

    def flat(self) -> list[_CubeNode[S]]:
        # List of nodes from start to head
        path = self
        nodes = [path.head]
        while path.prev is not None:
            path = path.prev
            nodes.append(path.head)
        return list(reversed(nodes))

    def support(self) -> set[CubeType]:
        # Set of cubes
        path = self
        cubes = {path.head.cube}
        while path.prev is not None:
            path = path.prev
            cubes.add(path.head.cube)
        return cubes

    def tail_support(self) -> set[CubeType]:
        # Set of cubes without head
        cubes = self.support()
        cubes.remove(self.head.cube)
        return cubes

    def future(self) -> tuple[frozenset[CubeType], _CubeNode[S] | None]:
        # It determines continuation of this path
        return frozenset(self.tail_support()), self.head

    def future_hash(self) -> int:
        # Hash of path future
        return hash(self.future())

    def complement_future_hash(self, all_cubes: Set[CubeType]):
        # future_hash of path that can be linked to self (other.head=self.head)
        complement_cubeset = all_cubes - self.support()
        complement_future = (frozenset(complement_cubeset), self.head)
        return hash(complement_future)

    def __str__(self) -> str:
        return ' -> '.join(str(node) for node in self.flat())


class CubePathTree[S]:
    # CubePathTree - knows how to grow tree of cube nodes using given rules.
    # Node is a pair (cube, state); cube is a usual dim-div-cube,
    # e.g. tuple with dim coords, in {0,..,div-1}
    # State is any hashable object. Grow semantics is defined by next_func.

    type CubeDeltaType = tuple[Literal[-1, 1, 0]]

    def __init__(self, dim: int, div: int,
            next_func: Callable[[CubeType, S], Iterable[tuple[CubeDeltaType, S]]],  # TODO: declare NextFuncType
            prev_func: Callable[[CubeType, S], Iterable[tuple[CubeDeltaType, S]]] | None = None):
        # next_func: sets the rule (cube, state) -> [(delta, new_state), ...] for growing paths
        # prev_func: used to prepare finish paths (optimization), grow backwards
        self.dim = dim
        self.div = div
        self.next_func = next_func
        self.prev_func = prev_func
        self._all_cubes = frozenset(itertools.product(range(div), repeat=dim))

    @classmethod
    def _init_path(cls, cube: CubeType, state: S) -> _CubePath:
        # Create path with initial node
        head = _CubeNode(cube=cube, state=state)
        return _CubePath(head=head, prev=None, length=1)

    def _continue_path(self, path: _CubePath[S], reverse: bool = False, cubeset: set[CubeType] | None = None) -> Iterable[_CubePath[S]]:
        # grow path one more node
        # reverse -  use prev_func to grow tree backwards
        # cubeset -  support of path (optimization!)
        N = self.div
        if cubeset is None:
            cubeset = path.support()

        cube = path.head.cube
        continue_func = self.prev_func if reverse else self.next_func
        for cube_delta, new_state in continue_func(path.head.cube, path.head.state):  # type: ignore
            new_cube = tuple(cj + dj for cj, dj in zip(cube, cube_delta))
            if any(nj < 0 or nj >= N for nj in new_cube) or (new_cube in cubeset):
                continue
            head = _CubeNode(cube=new_cube, state=new_state)
            yield _CubePath(head, length=path.length + 1, prev=path)

    def _depth_search(self, paths: Iterable[_CubePath[S]], length: int) -> Iterable[_CubePath[S]]:
        # Depth-first search for paths.
        # paths   -  starting paths
        # length  -  yield paths of that length
        todo = [(path, path.support()) for path in paths]
        while todo:
            path, cubeset = todo.pop()
            if path.length == length - 1:  # almost done
                yield from self._continue_path(path, cubeset=cubeset)
            elif path.length < length - 1:
                for new_path in self._continue_path(path, cubeset=cubeset):
                    new_cubeset = cubeset.copy()
                    new_cubeset.add(new_path.head.cube)
                    todo.append((new_path, new_cubeset))

    def _width_search(self, paths: list[_CubePath[S]], max_steps: int, max_count: int | None = None, reverse: bool = False) -> list[_CubePath[S]]:
        # returns paths from one step to maintain equal length
        curr_paths = paths
        for i in range(max_steps):
            new_paths: list[_CubePath[S]] = []
            for path in curr_paths:
                new_paths += self._continue_path(path, reverse=reverse)
                if max_count is not None and len(new_paths) > max_count:
                    logger.debug('_width_search stopped due to max_count %d', max_count)
                    return curr_paths
            curr_paths = new_paths
            logger.debug('_width_search: step %d, expanded to %d', i+1, len(curr_paths))
        return curr_paths

    def grow(self,
            start: Iterable[tuple[CubeType, S]],
            finish: Iterable[tuple[CubeType, S]],
            start_max_count: int = 100,
            finish_max_count: int | None = 10**6
        ) -> Iterable[list[tuple[CubeType, S]]]:
        # grow cube paths that fill div-grid
        # start - pairs (cube, state) for first node in path
        # finish - pairs (cube, state) for last node in path
        # start_max_count -- number of width-search-expanded start paths
        # finish_max_count -- number of width-search-expanded finish paths (maximize to speed up search)
        # yields flat paths: lists of pairs (cube, state) from start to end

        N = self.div
        d = self.dim
        max_steps = (N**d // 2) - 1  # force that start and finish paths do not intersect

        start_paths = [self._init_path(cube, state) for cube, state in start]
        start_paths = self._width_search(start_paths, max_steps=max_steps, max_count=start_max_count)
        if not start_paths:
            logger.debug('start: cannot find start paths')
            return

        finish_paths = [self._init_path(cube, state) for cube, state in finish]
        if finish_max_count is not None:
            if self.prev_func is None:
                raise Exception("prev_func not defined, can't expand finish")
            finish_paths = self._width_search(finish_paths, max_steps=max_steps, max_count=finish_max_count, reverse=True)
        if not finish_paths:
            logger.debug('finish: cannot find paths')
            return
        finish_by_hash = defaultdict(list)
        for path in finish_paths:
            finish_complemented_hash = path.complement_future_hash(self._all_cubes)
            finish_by_hash[finish_complemented_hash].append(path)
        mid_length = N**d - finish_paths[0].length + 1  # +1 as mid & finish = head

        start_path_groups = list(self._group_paths_by_future(start_paths))
        logger.debug(
            'depth_search width %d: start %d paths, %d groups => finish %d paths, %d hashes',
            mid_length - start_paths[0].length, len(start_paths), len(start_path_groups), len(finish_paths), len(finish_by_hash)
        )

        found_count = 0
        for cnt, paths in enumerate(start_path_groups):
            logger.debug('processing start: %d of approx %d, found: %d', cnt + 1, len(start_path_groups), found_count)
            for mid_path in self._depth_search([paths[0]], length=mid_length):
                fin_paths = finish_by_hash.get(mid_path.future_hash(), [])
                fin_paths = [path for path in fin_paths if self._check_glue(mid_path, path)]
                for fin_path in fin_paths:
                    for start_path in paths:
                        yield self._glue_paths(start_path, mid_path, fin_path)
                        found_count += 1

    @staticmethod
    def _group_paths_by_future(paths: Iterable[_CubePath[S]]) -> Iterable[list[_CubePath[S]]]:
        # first we group by future_hash (optimization)
        pdict = defaultdict(list)
        for path in paths:
            pdict[path.future_hash()].append(path)
        for paths_with_hash in pdict.values():
            future2paths = defaultdict(list)
            for path in paths_with_hash:
                future2paths[path.future()].append(path)
            yield from future2paths.values()

    @staticmethod
    def _check_glue(beg: _CubePath[S], end: _CubePath[S]) -> bool:
        # Check that end path is a continuation
        return beg.head == end.head and (not end.tail_support() & beg.tail_support())

    @staticmethod
    def _glue_paths(start: _CubePath[S], mid: _CubePath[S], finish: _CubePath[S]) -> list[tuple[CubeType, S]]:
        # Glue mid -- continuation of start_future(!) -- to fin path
        finish_reversed = list(reversed(finish.flat()))
        pdata = start.flat() + mid.flat()[start.length:] + finish_reversed[1:]
        return [(node.cube, node.state) for node in pdata]
