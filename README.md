# Description

Python package "peano" implements regular peano fractal curves on d-dimensional cube.

Required python packages:
- sympy (Rational module)
- python-sat

Tested on python 3.5 and 3.8.

This repo is mostly frozen, only fixes will appear.


# User and developer guide

Entry point is the script search.py: try it!

Some examples of Peano curves are in test/examples.py

Operator usage:
b * X means action of b (usually, base_map) on object X (curve, path, etc)
~X means inverted object: "inverted map" or "reversed curve"

We prefer to use generators if possible.

We use immutable hashable objects if possible (e.g., tuples intead of lists).

Common names:
* dim -- dimension of the curve image: [0,1]->[0,1]^d
* div -- number of divisions for each coordinate
* pcount -- number of patterns (monofractal: 1, bifractal: 2, etc)
* cube -- d-tuple of ints, all coordinates range from 0 to div-1
* cnum -- index of cube, ranges from 0 to genus-1
* pnum -- index of pattern, ranges from 0 to pcount-1
* genus = div**dim -- genus of the curve
