import itertools
import logging
from collections import namedtuple, defaultdict


class CubeNode(namedtuple('CubeNode', ['cube', 'state'])):
    """
    Node of the partial path in the cube lattice.

    cube  --  cube position
    state --  state (i.e., position inside the cube), any hashable object
    """
    def __str__(self):
        return '{} @{}'.format(self.cube, self.state)


class NodePath(namedtuple('NodePath', ['head', 'len', 'prev'])):
    """
    Path of nodes, with usual head+link_to_tail representation.

    head  --  head node (CubeNode)
    len   --  path length
    prev  --  link to tail path
    """

    @classmethod
    def get_empty(cls):
        """Empty path."""
        return cls(head=None, len=0, prev=None)

    def flat(self):
        """List of nodes from start to head."""
        nodes = []
        path = self
        while path.head is not None:
            nodes.append(path.head)
            path = path.prev
        return list(reversed(nodes))

    def __invert__(self):
        """Reversed path, rather inefficient."""
        curr = self.get_empty()
        for node in reversed(self.flat()):
            new_path = type(self)(head=node, prev=curr, len=curr.len + 1)
            curr = new_path
        return curr

    def get_cubes(self, skip_head=False):
        length = self.len
        path = self
        if skip_head:
            length -= 1
            path = path.prev
        cubes = [None] * length
        for i in range(length):
            cubes[i] = path.head.cube
            path = path.prev
        return cubes

    def support(self):
        """Set of cubes."""
        return set(self.get_cubes())

    def tail_support(self):
        """Set of cubes without head"""
        return set(self.get_cubes(skip_head=True))

    def future(self):
        """It determines continuation of this path."""
        return frozenset(self.tail_support()), self.head

    def phash(self):
        """Hash of path future."""
        return hash(self.future())

    def __str__(self):
        return ' -> '.join([str(node) for node in self.flat()])


class NodePathTree:
    """
    NodePathTree instance knows how to grow tree of nodes using given rules.
    Each node's state is just some object, all grow semantics is in the next_func.
    """

    def __init__(self, dim, div, next_func, prev_func=None):
        """
        next_func: sets the rule for growing:  head -> list of pairs (cube_delta, new_state)
        prev_func: used to prepare finish paths (optimization), grow backwards
        """
        self.dim = dim
        self.div = div
        self.next_func = next_func
        self.prev_func = prev_func
        self.all_cubes = frozenset(itertools.product(range(div), repeat=dim))

    @classmethod
    def init_path(cls, cube, state):
        """Create path with initial node."""
        head = CubeNode(cube=cube, state=state)
        return NodePath(head=head, prev=NodePath.get_empty(), len=1)

    def get_complement_phash(self, path):
        """Phash of any path that can be linked to given path."""
        complement_cubeset = self.all_cubes - path.support()
        complement_future = (frozenset(complement_cubeset), path.head)
        return hash(complement_future)

    def continue_path(self, path, reverse=False, cubeset=None):
        """
        Args:
            path    -  subj
            reverse -  use prev_func to grow tree backwards
            cubeset -  support of path (optimization!)
        """
        N = self.div
        if cubeset is None:
            cubeset = path.support()

        cube = path.head.cube
        continue_func = self.prev_func if reverse else self.next_func
        for cube_delta, new_state in continue_func(path.head):
            new_cube = tuple(cj + dj for cj, dj in zip(cube, cube_delta))
            if any(nj < 0 or nj == N for nj in new_cube) or new_cube in cubeset:
                continue
            head = CubeNode(cube=new_cube, state=new_state)
            yield NodePath(head, len=path.len + 1, prev=path)

    def depth_search(self, paths, max_len, finish_phashes=None):
        """Depth-first search for paths.
        Params:
            paths           -  starting paths
            max_len         -  subj
            finish_phashes  -  finish path hashes (set or dict), path should intersect finish by head
        """
        todo = [(path, path.support()) for path in paths]
        while todo:
            path, cubeset = todo.pop()
            if path.len == max_len - 1:  # almost done
                for path in self.continue_path(path, cubeset=cubeset):
                    if finish_phashes is not None and path.phash() not in finish_phashes:
                        continue
                    yield path
            else:
                for np in self.continue_path(path, cubeset=cubeset):
                    new_cubeset = cubeset.copy()
                    new_cubeset.add(np.head.cube)
                    todo.append((np, new_cubeset))

    def width_search(self, paths, max_steps, max_count=None, reverse=False):
        """Width-search for given starting paths."""
        curr_paths = paths
        stop = False
        for i in range(max_steps):
            new_paths = []
            for path in curr_paths:
                new_paths += self.continue_path(path, reverse=reverse)
                if max_count is not None and len(new_paths) > max_count:
                    stop = True
                    break
            if stop:
                break
            curr_paths = new_paths
            logging.debug('width_search: step %d, expanded to %d' % (i+1, len(curr_paths)))
        return curr_paths

    @classmethod
    def group_paths(cls, paths):
        pdict = defaultdict(list)
        # group start paths by future because it determines continuation
        # we use that phash=hash(future), key=future was too expensive
        for path in paths:
            pdict[path.phash()].append(path)

        for phash_paths in pdict.values():
            future2paths = defaultdict(list)
            for path in phash_paths:
                future2paths[path.future()].append(path)
            yield from future2paths.values()

    def grow(self, paths, ends=None, start_max_count=100, finish_max_count=10**6):
        """
        paths           -- start paths to grow
        ends            -- path ends (if None, end anywhere), i.e. result path.head in end heads set
        start_max_count -- number of start paths
        finish_max_count -- number of finish paths (maximize given RAM)
        """
        N = self.div
        d = self.dim
        max_steps = (N**d // 2) - 1  # force that start and finish paths do not intersect

        start_paths = self.width_search(paths, max_steps=max_steps, max_count=start_max_count)
        if not start_paths:
            logging.debug('start: cannot find start paths')
            return
        logging.debug('start: width: %d, paths: %d', start_paths[0].len, len(start_paths))

        if ends:
            if any(path.len != ends[0].len for path in ends):
                raise Exception("All end paths must have same length (optimization)!")

            finish_paths = [~path for path in ends]
            if finish_max_count is not None:
                finish_paths = self.width_search(finish_paths, max_steps=max_steps, max_count=finish_max_count, reverse=True)
            if not finish_paths:
                logging.debug('finish: cannot find paths')
                return
            finish_width = finish_paths[0].len

            logging.debug('finish: width %d, paths: %d', finish_width, len(finish_paths))
            finish = defaultdict(list)
            for path in finish_paths:
                phash = self.get_complement_phash(path)
                finish[phash].append(path)
            depth_kwargs = {'max_len': N**d - finish_width + 1, 'finish_phashes': finish}
        else:
            depth_kwargs = {'max_len': N**d}

        found_count = 0
        for cnt, paths in enumerate(self.group_paths(start_paths)):
            if (cnt + 1) % 100 == 0:
                logging.debug('processing start: %d of approx %d, found: %d', cnt + 1, len(start_paths), found_count)
            for mid_path in self.depth_search([paths[0]], **depth_kwargs):
                if ends is not None:
                    fin_paths = finish[mid_path.phash()]
                else:
                    fin_paths = [NodePath.get_empty()]
                for fin_path in fin_paths:
                    if not self.check_glue(mid_path, fin_path):  # start and mid are correctly glued
                        continue
                    for start_path in paths:
                        path_data = self.glue_paths(start_path, mid_path, fin_path)
                        if path_data is None:
                            continue
                        yield path_data
                        found_count += 1

    @staticmethod
    def check_glue(beg, end):
        """Check that end path is a continuation."""
        if end.head is None:  # empty path
            return True
        if beg.head != end.head:
            return False
        if end.tail_support() & beg.tail_support():
            return False
        return True

    @staticmethod
    def glue_paths(start, mid, fin):
        """Glue mid path (continuation of start) to fin path"""
        rev_fin = list(reversed(fin.flat()))
        return start.flat() + mid.flat()[start.len:] + rev_fin[1:]
