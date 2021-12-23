# Description

Python package "peano" implements regular peano fractal curves on d-dimensional cube.
This code is used to obtain results described in the paper:
https://arxiv.org/abs/2103.07344

Computations from the paper may be reproduced using the [a relative link](search.py) script:

Computation 1 (Spring 3d curve):
* $ ./search.py --dim 3 --pcount 2 --div 2 --facet-gated --metric l2_squared --rel-tol-inv 1000000
* time: ~4 min
Note that output curves differ from the description in the paper.
See also ``get_spring_curve'' method in [a relative link](tests/examples.py).

Computation 2 (4d curve):
* $ ./search.py --dim 4 --pcount 1 --div 2 --facet-gated --metric l2 --rel-tol-inv 100000
* time: ~45 min

* YE curve search
* $ ./search.py --dim 2 --pcount 1 --div 5 --gates '(0,0)->(0,1)' --max-cdist 1 --metric l2 --rel-tol-inv 1000000
* time: ~1.5min

See also [a relative link](YE-proof.py) for a full computer-assisted proof of YE curve minimality.

In our test examples we use data from several papers:

Haverkort, Walderveen, ``Locality and bounding-box quality of two-dimensional
space-filling curves'', 2010

Haverkort, ``An inventory of three-dimensional Hilbert space-filling curves'', 2011
https://arxiv.org/abs/1109.2323

Korneev, Shchepin, ``L-infty-Locality of Three-Dimensional Peano Curves'', 2019

# User guide

Requirements: see requirements.txt

Tested on python 3.5 and 3.8.

This repo is mostly frozen, only fixes will appear.

Some other examples of search.py:
* facet-gates search:
    * $ ./search.py --dim 3 --pcount 2 --div 2 --facet-gated --output-gates
    * time: ~20s, result: 35 gates
* bifractal gates search:
    * $ ./search.py --dim 2 --pcount 2 --div 2 --output-gates
    * time: ~2m30s, result: 9 gates

Some Peano curves may be found in: tests/examples.py.

# Developer guide

Operator usage:
* b * X means action of b (usually, base_map) on object X (curve, path, etc)
* b\*\*(-1) is group inverse of b, i.e. c such that b\*c=c\*b=id
* ~X means time reversal (in base map or curve)

Common names and notations:
* dim -- dimension of the curve image: [0,1]->[0,1]^d
* div -- number of divisions for each coordinate
* pcount -- number of patterns (monofractal: 1, bifractal: 2, etc)
* cube -- d-tuple of ints, all coordinates range from 0 to div-1
* cnum -- index of cube, ranges from 0 to genus-1
* pnum -- index of pattern, ranges from 0 to pcount-1
* genus = div\*\*dim -- genus of the curve
* proto -- prototype (sequence of cubes)
* path -- prototype with entrances/exits in each cube

We prefer to use generators if possible.

We use immutable hashable objects if possible (e.g., tuples intead of lists).
