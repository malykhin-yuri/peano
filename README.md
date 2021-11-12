# Description

Python package "peano" implements regular peano fractal curves on d-dimensional cube.
This code is used to obtain results described in the paper:
https://arxiv.org/abs/2103.07344

Requirements: see requirements.txt

Tested on python 3.5 and 3.8.

This repo is mostly frozen, only fixes will appear.

# User and developer guide

Entry point is the script search.py: try it!
Results stated in the arXiv paper:
* YE curve:
    * $ ./search.py --dim 2 --pcount 1 --div 5 --gates '(0,0)->(0,1)' --max-cdist 1 --metric l2 --rel-tol-inv 1000000 --output-curves
    * time: ~3m30s, result: ye curve
* Spring curve:
    * $ ./search.py --dim 3 --pcount 2 --div 2 --facet-gated --metric l2_squared --output-curves
    * time: ~7min
* 4D facet-gated curve:
    * $ ./search.py --dim 4 --pcount 1 --div 2 --facet-gated --metric l2 --output-curves
    * time: ?

Some other examples:
* facet-gates search:
    * $ ./search.py --dim 3 --pcount 2 --div 2 --facet-gated --output-gates
    * time: ~20s, result: 35 gates
* bifractal gates search:
    * $ ./search.py --dim 2 --pcount 2 --div 2 --output-gates
    * time: ~2m30s, result: 9 gates

Some Peano curves may be found in: tests/examples.py.

See also YE-proof.py for a computer-assisted proof of YE curve minimality.

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
* genus = div**dim -- genus of the curve
* proto -- prototype (sequence of cubes)
* path -- prototype with entrances/exits in each cube

We prefer to use generators if possible.

We use immutable hashable objects if possible (e.g., tuples intead of lists).
