# Description

Python package "peano" implements regular peano fractal curves on d-dimensional cube
and the algorithm for the search of minimal dilation curves.

We use modern typing syntax that requires Python 3.12+
(use the [da691](https://github.com/malykhin-yuri/peano/tree/da691a2c2f262f559445be2d2ede70563c3ea1f1) version for Python 3.8).
The dependencies are listed here: [requirements](requirements.txt).
We rely on `python-sat` module for SAT-solvers and `quicktions` module for fast `Fractions`.

# Published computations

This code was used to obtain results described in the paper
"Search of fractal space-filling curves with minimal dilation", see [paper-search](paper-search.md).

# Tutorial

## Peano multifractals

Peano curve is a surjective map [0,1]->[0,1]^d. We consider the class of Peano
multifractals which is precisely defined in the [paper](paper-search.md). Divide [0,1]^d to s^d
cubes with side 1/s; let us call them s-cubes.
We say that a tuple of Peano curves (p_1,...,p_m) is 
a **Peano multifractal** if each curve p_j on each s-cube is similar
to one of the curves {p_i}. Similarity is defined by a **base_map** which
is an isometry of [0,1]^d and isometry of [0,1] (time reversal or identity map).

So, multifractal curve has three main characteristics: **dimension** (i.e., d),
**multiplicity** (i.e., m) -- number of patterns (m=1 is called a monofractal, m=2 a bifractal, etc)
and **genus** -- number of first-order fractions; genus equals d^s.

To completely define a multifractal curve, we must specify m **patterns**. For each
pattern two things are required. First is a **prototype** -- a sequence of
sub-cubes of [0,1]^d that defines curve order. Next, for each cube in the
prototype we define **spec** -- i.e. we specify what base_map and what
pattern may be used to obtain sub-curve on that cube.

Fairly complete collection of the curves and some information about them
is presented in our **Zoo**: [zoo.py](peano/zoo.py).
There we use compact notation to define multifractals; it is briefly exaplained
in `Proto.parse` method in [paths](peano/paths.py) module and `Spec.parse`,
`BaseMap.parse` methods in [base_maps](peano/base_maps.py) module.

In these examples we use data from several papers [HW10], [H11], [KS18]
(see the "Bibliography" section below).

Here are some common notations used in code:
* `dim` -- dimension of the curve (also called `d`)
* `div` -- number of divisions for each coordinate (called `s` above)
* `mult` -- multiplicity of the curve, i.e. number of patterns
* `cube` -- div-cube {c_i/div <= x_i <= (c_i+1)/div}, coded as d-tuple of integers (c_0,...,c_{d-1}), where 0 <= c_i < div
* `cnum` -- index of cube, ranges from 0 to genus-1
* `pnum` -- index of pattern, ranges from 0 to mult-1
* `genus` = div\*\*dim -- genus of the curve
* `proto` -- prototype (sequence of cubes)
* `path` -- prototype with entrances/exits in each cube
* b * X means action of b (usually, base_map) on object X (curve, path, etc)
* b\*\*(-1) is group inverse of b, i.e. c such that b\*c=c\*b=id
* ~X means time reversal (in base map or curve)

Curves are implemented in the [curves](peano/curves.py) module. Start with the
class Curve there. There are some basic operations that one can want to do with
curves:
* `get_entrance`, `get_exit`
* `get_vertex_moments` -- get moments of time when curve visits cube vertices 
* apply isometry to the curve: b * C, where b is an ``BaseMap`` instance

Specs and base maps are defined in the corresponding module [base_maps](peano/base_maps.py).

See also examples in [tests](tests). 

## Dilation

The main purpose of this module is to search curves with minimal dilation. It is
defined as

```max_{0<=s<t<=1} ||p(s)-p(t)||^d / |s-t|```

Here ||.|| is some norm, e.g. l_p-norm.
This maximum is always attained. To work with the dilation you need
some theory that is given in our cited paper and other papers of E.V.Shchepin.
E.g. you have to use `curves.Junctions`. The main tool to
estimate dilation and search of minimal-dilation curves is the class `Estimator`
from the [dilation](peano/dilation.py) module. See also its documentation.

Estimator uses SAT-solvers to estimate the minimum dilation of the class of all
curves with given pointed prototype (prototype with fixed entrances/exits in
each fraction). That's why we need the class `FuzzyCurve`. The algorithm is
described in our paper. All possible pointed prototypes are generated in two
steps. First, using [gates](peano/gates.py) module we
generate all possible **gates** (entrance/exit configurations). Second, using
[paths](peano/paths.py) we generate all possible pointed prototypes (in our
code, we use alias `path` for it) with given gates.

All described machinery is combined in the [search](search.py) script. We already
mentioned a couple of cases, here are some other:
* facet-gates search:
    * $ ./search.py --dim 3 --mult 2 --div 2 --facet-gated --output-gates
    * time: ~20s, result: 35 gates
* bifractal gates search:
    * $ ./search.py --dim 2 --mult 2 --div 2 --output-gates
    * time: ~2m30s, result: 9 gates

## Development

The library is open source and you are encouraged to fork and extend it.
Note that we prefer to use generators if possible and we use immutable
hashable objects if possible (e.g., tuples intead of lists).

Please note the essential restriction of the `Estimator`: fuzzy curves must
have independent fractions, i.e. all possible combinations of variants for each
fraction must give a correct curve. This may be avoided somehow, but it will
require changes in the SAT-adapter.

Useful:
* static type check:
    $ python -m mypy peano tests  # ensure using venv python
* unit tests:
    $ python -m unittest discover -f

# Additional links
Another python library about Peano curves, with useful table of records: https://github.com/Anton495/space-filling-curves

# Bibliography

[H11] H. Haverkort, An inventory of three-dimensional Hilbert space-ﬁlling curves: E-print, 2011. arXiv: 1109.2323v2

[HW10] H. Haverkort, F. Walderveen, "Locality and bounding-box quality of two-dimensional
space-ﬁlling curves", *Comp. Geom. Th. Appl.*, **43** (2010), 131--147.

[KS18] A.A. Korneev, E.V. Shchepin, *Proceedings of the Steklov Institute of Mathematics*, **302** (2018), 217--249.

[MS23] Yu. Malykhin, E. Shchepin,
"Search of fractal space-filling curves with minimal dilation",
*Discrete Comput. Geom.*, **70** (2023), 189--213.

[SB08] E.V. Shchepin, K.E. Bauman, "Minimal Peano Curve", Proc. Steklov Inst. Math., **263** (2008), 236--256.

[Tok10] С.С. Токарев. Дипломная работа. 2010.
