import sys
import unittest
import logging
import itertools

from peano._cube_path_trees import CubePathTree


class Next:
    def __init__(self, dim, hyper=False):
        deltas = list(itertools.product((-1, 0, 1), repeat=dim))
        if hyper:
            deltas = [dt for dt in deltas if sum(abs(x) for x in dt) == 1]
        self.deltas = deltas

    def __call__(self, cube, state):
        return [(delta, None) for delta in self.deltas]


class TestCubePathTree(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.INFO, stream=sys.stdout)

    def test_king(self):
        # https://oeis.org/A272445 - Numbers of paths for moving a king from a corner to the opposite one
        # in a n X n chessboard, provided that each cell must be reached exactly once
        true_results = [0, 1, 2, 30, 4942]
        dim = 2
        for N in range(2, 4):
            tree = CubePathTree(dim=dim, div=N, next_func=Next(dim))
            start_cube = (0, 0)
            finish_cube = (N-1, N-1)
            found_paths = list(tree.grow(start=[(start_cube, None)], finish=[(finish_cube, None)], finish_max_count=None))
            self.assertEqual(len(found_paths), true_results[N])

        for N in range(4, 5):
            tree = CubePathTree(dim=dim, div=N, next_func=Next(dim), prev_func=Next(dim))
            start_cube = (0, 0)
            finish_cube = (N-1, N-1)
            found_paths = list(tree.grow(start=[(start_cube, None)], finish=[(finish_cube, None)], finish_max_count=200))
            self.assertEqual(len(found_paths), true_results[N])

    def test_grow_proto(self):
        state = None
        dim = 2
        div = 3
        start = (0,) * dim, state
        end = (div-1,) * dim, state
        ptree = CubePathTree(dim=dim, div=div, next_func=Next(dim, hyper=True))
        paths3 = list(ptree.grow([start], [end], finish_max_count=None))
        assert len(paths3) == 2

        ptree2 = CubePathTree(dim=dim, div=div, next_func=Next(dim, hyper=True), prev_func=Next(dim, hyper=True))
        paths4 = list(ptree2.grow([start], [end], finish_max_count=100))
        assert len(paths4) == 2
