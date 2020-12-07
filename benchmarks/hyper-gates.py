#!/usr/bin/env python3
"""
Benchmark: find possible gates for 2d bifractals of genus 4.

Time: 18s
Result: 35 good gates
"""

import logging
import sys
sys.path.append('.')

from peano.gate_utils import GatesGenerator

logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
    list(GatesGenerator(dim=3, div=2, pattern_count=2, hyper=True).gen_gates())
