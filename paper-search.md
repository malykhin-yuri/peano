# Paper "Search of fractal space-filling curves with minimal dilation"

Preprint: https://arxiv.org/abs/2103.07344

Computations from the paper may be reproduced using the [search.py](search.py) script:

Computation 1 (Spring 3d curve):
* $ ./search.py --dim 3 --mult 2 --div 2 --facet-gated --metric l2_squared --rel-tol-inv 1000000
* time: ~4 min

Note that output curves may differ from the description in the paper (by some base\_map).
See also method `get_spring_curve` in [examples](tests/examples.py).

Computation 2 (4d curve):
* $ ./search.py --dim 4 --mult 1 --div 2 --facet-gated --metric l2 --rel-tol-inv 100000
* time: ~45 min

See an example of an output curve in the method `get_4d_facet_gated_curve` in examples.

YE curve search:
* $ ./search.py --dim 2 --mult 1 --div 5 --gates '(0,0)->(0,1)' --max-cdist 1 --metric l2 --rel-tol-inv 1000000
* time: ~1.5min

See also [YE-proof](YE-proof.py) for a full computer-assisted proof of YE curve minimality.

In the [logs](logs) directory we put the output of three mentioned computations.

