# Description

Python package "peano" implements regular peano fractal curves on d-dimensional cube
and the algorithm for the search of minimal dilation curves.

Code was tested on python 3.8 on Ubuntu Linux, see also the [requirements](requirements.txt).

This repo is mostly frozen, only fixes will appear.

# Published computations

This code is used to obtain results described in the paper:
https://arxiv.org/abs/2103.07344

Computations from the paper may be reproduced using the [search.py](search.py) script:

Computation 1 (Spring 3d curve):
* $ ./search.py --dim 3 --pcount 2 --div 2 --facet-gated --metric l2_squared --rel-tol-inv 1000000
* time: ~4 min

Note that output curves may differ from the description in the paper (by some base\_map).
See also method ``get_spring_curve'' in [examples](tests/examples.py).

Computation 2 (4d curve):
* $ ./search.py --dim 4 --pcount 1 --div 2 --facet-gated --metric l2 --rel-tol-inv 100000
* time: ~45 min

See an example of an output curve in the method ``get_4d_facet_gated_curve'' in [examples](tests/examples.py).

* YE curve search
* $ ./search.py --dim 2 --pcount 1 --div 5 --gates '(0,0)->(0,1)' --max-cdist 1 --metric l2 --rel-tol-inv 1000000
* time: ~1.5min

See also [YE-proof](YE-proof.py) for a full computer-assisted proof of YE curve minimality.

In the [logs](logs) directory we put the output of three mentioned computations.

In our test examples we use data from several papers:

Haverkort, Walderveen, ``Locality and bounding-box quality of two-dimensional
space-filling curves'', 2010

Haverkort, ``An inventory of three-dimensional Hilbert space-filling curves'', 2011
https://arxiv.org/abs/1109.2323

Korneev, Shchepin, ``L-infty-Locality of Three-Dimensional Peano Curves'', 2019

# Tutorial

## Peano multifractals

Peano curve is a surjective map [0,1]->[0,1]^{dim}. We consider the class of Peano
multifractals which is precisely defined in our paper. Briefly speaking, each
multifractal is a tuple of Peano curves (p_1,...,p_{mult}), each curve p_j is divided on cubic
fractions -- we divide [0,1]^d on sub-cubes with side 1/{div}, where div is an
integer; on each fraction p_j is isometric to one of the curves {p_i}.

Multifractal curve has three main characteristics: **dimension** (i.e., dim),
**multiplicity** -- number of patterns (1 for monofractal, 2 for bifractal, etc)
and **genus** -- number of first-order fractions. Genus equals dim^div.

To completely define a multifractal curve, we must specify two things. First,
**prototypes** for each pattern. A prototype is a sequence of sub-cubes of
[0,1]^d that defines curve order. Next, for each cube in each prototype we
define **spec** -- i.e. we specify what **base_map** and what pattern may be
used to obtain sub-curve on that cube. A base map is a pair of cube isometry and
[0,1] isometry (time reversal or identity map).

Many examples of the curves may be found in tests: [examples.py](tests/examples.py).
There we use compact notation to define multifractals. It is briefly examplained
in `Proto.parse` method in [paths](peano/paths.py) module; `Spec.parse` and
`BaseMap.parse` methods in [base_maps](peano/base_maps.py) module.

Here are some common notations used in code:
* dim -- dimension of the curve
* div -- number of divisions for each coordinate
* mult -- multiplicity of the curve, i.e. number of patterns
* cube -- d-tuple of integers, all coordinates range from 0 to div-1
* cnum -- index of cube, ranges from 0 to genus-1
* pnum -- index of pattern, ranges from 0 to mult-1
* genus = div\*\*dim -- genus of the curve
* proto -- prototype (sequence of cubes)
* path -- prototype with entrances/exits in each cube
* b * X means action of b (usually, base_map) on object X (curve, path, etc)
* b\*\*(-1) is group inverse of b, i.e. c such that b\*c=c\*b=id
* ~X means time reversal (in base map or curve)

Curves are implemented in the [curves](peano/curves.py) module. Start with the
class Curve there. There are some basic operations that one can want to do with
curves:
* `get_entrance`, `get_exit`
* `get_vertex_moments` -- get moments of time when curve visits cube vertices 
* apply isometry to the curve: b * C, where b is an ``BaseMap`` instance

See also examples in tests.

## Dilation

The main purpose of this module is to search curves with minimal dilation. It is
defined as
```max ||p(s)-p(t)||^{dim} / |s-t|```
This maximum (over s,t in [0,1]) is always attained. To estimate it you require
some theory that is given in our cited paper and other papers of E.V.Shchepin.
You have to use `Junctions` from [curves](peano/curves.py). The main tool to
estimate dilation and search of minimal-dilation curves is the class `Estimator`
from [dilation](peano/dilation.py) module. See also its documentation.

Estimator uses SAT-solvers to estimate the minimum dilation of the class of all
curves with given pointed prototype (prototype with fixed entrances/exits in
each fraction). The algorithm is described in our paper. All possible pointed prototypes
are generated in two steps. First, using [gates](peano/gates.py) module we
generate all possible **gates** (entrance/exit configurations). Second, using
[paths](peano/paths.py) we generate all possible pointed prototypes (in our
code, we use alias `path` for it) with given gates.

All described logic is combined in the [search](search.py) script. We already
mentioned a couple of cases, here are some other:
* facet-gates search:
    * $ ./search.py --dim 3 --pcount 2 --div 2 --facet-gated --output-gates
    * time: ~20s, result: 35 gates
* bifractal gates search:
    * $ ./search.py --dim 2 --pcount 2 --div 2 --output-gates
    * time: ~2m30s, result: 9 gates

## Development
The library is open source and you are encouraged to extend it.
Note that we prefer to use generators if possible and we use immutable
hashable objects if possible (e.g., tuples intead of lists).
