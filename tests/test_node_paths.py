import unittest
import itertools

from peano.node_paths import NodePathTree


class Next:
    def __init__(self, dim, hyper=False):
        deltas = list(itertools.product((-1, 0, 1), repeat=dim))
        if hyper:
            deltas = [dt for dt in deltas if sum(abs(x) for x in dt) == 1]
        self.deltas = deltas
        self.state = None

    def __call__(self, _):
        return [(delta, self.state) for delta in self.deltas]


class TestNodePathTree(unittest.TestCase):
    def test_grow_proto(self):
        state = None
        dim = 2
        start = NodePathTree.init_path((0,) * dim, state)

        paths = list(NodePathTree(dim=dim, div=2, next_func=Next(dim)).grow([start]))
        assert len(paths) == 6

        paths2 = list(NodePathTree(dim=dim, div=2, next_func=Next(dim, hyper=True)).grow([start]))
        assert len(paths2) == 2

        dim = 2
        div = 3
        start = NodePathTree.init_path((0,) * dim, state)
        end = NodePathTree.init_path((div-1,) * dim, state)
        ptree = NodePathTree(dim=dim, div=div, next_func=Next(dim, hyper=True))
        paths3 = list(ptree.grow([start], ends=[end], finish_max_count=None))
        assert len(paths3) == 2

        ptree2 = NodePathTree(dim=dim, div=div, next_func=Next(dim, hyper=True), prev_func=Next(dim, hyper=True))
        paths4 = list(ptree2.grow([start], ends=[end], finish_max_count=100))
        assert len(paths4) == 2
