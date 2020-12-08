import unittest

from peano.gate_utils import GatesGenerator


class TestGates(unittest.TestCase):
    def test_gen_gates(self):
        assert len(list(GatesGenerator(dim=3, div=2, pcount=1).gen_gates())) == 6

    def test_gen_gates(self):
        assert len(list(GatesGenerator(dim=2, div=2, pcount=2, hyper=True).gen_gates())) == 1
        assert len(list(GatesGenerator(dim=2, div=3, pcount=2, hyper=True).gen_gates())) == 2
        assert len(list(GatesGenerator(dim=2, div=2, pcount=3, hyper=True).gen_gates())) == 6
        assert len(list(GatesGenerator(dim=3, div=2, pcount=2, hyper=True).gen_gates())) == 35
