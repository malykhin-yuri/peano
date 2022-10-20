import logging
import pprint

from peano.utils import ratio_l2
from peano.subsets import Link
from peano.base_maps import BaseMap
from peano.paths import PathsGenerator
from peano.curves import PathFuzzyCurve
from peano.dilation import Estimator
from tests.examples import get_ye_curve

from search import run_search


DIAG_LINK = Link.parse_gates('(0,0)->(1,1)')
SIDE_LINK = Link.parse_gates('(0,0)->(1,0)')
MEDIAN_LINK = Link.parse_gates('(0,0)->(1,1/2)')


def check_genus5_non_ye():
    YE_proto = get_ye_curve().proto
    YE_protos = set(bm * YE_proto for bm in BaseMap.gen_base_maps(dim=2))
    paths_gen = PathsGenerator(dim=2, div=5, links=(SIDE_LINK,), max_cdist=1)
    paths_list = list(paths_gen.generate_paths(std=True))
    print('got paths:', len(paths_list))
    paths_list = [paths for paths in paths_list if paths[0].proto not in YE_protos]
    print('got non-YE paths:', len(paths_list))

    estimator = Estimator(ratio_l2)

    result = estimator.estimate_dilation_sequence(
        [PathFuzzyCurve.init_from_paths(paths) for paths in paths_list],
        rel_tol_inv=200,
    )
    result.update(estimator.get_info())
    pprint.pprint(result)
    print('lower bound:', float(result['lo']))
    print('upper bound:', float(result['up']))


def check_genus5_ye_proto():
    YE_curve = get_ye_curve()
    YE_proto = YE_curve.proto
    YE_protos = set(bm * YE_proto for bm in BaseMap.gen_base_maps(dim=2))
    paths_gen = PathsGenerator(dim=2, div=5, links=(SIDE_LINK,), max_cdist=1)
    paths_list = list(paths_gen.generate_paths(std=True))
    print('got paths:', len(paths_list))
    paths_list = [paths for paths in paths_list if paths[0].proto in YE_protos]
    print('got YE paths:', len(paths_list))
    assert len(paths_list) == 1

    test_pcurves = []
    ye_pcurve = get_ye_curve().forget()
    for cnum, ye_spec in enumerate(YE_curve.specs):
        allowed = ye_pcurve.gen_allowed_specs(pnum=0, cnum=cnum)
        test_pcurves += [ye_pcurve.specify(pnum=0, cnum=cnum, spec=sp) for sp in allowed if sp != ye_spec]

    estimator = Estimator(ratio_l2)

    result = estimator.estimate_dilation_sequence(
        test_pcurves,
        rel_tol_inv=1000,
    )
    result.update(estimator.get_info())
    pprint.pprint(result)
    print('lower bound:', float(result['lo']))
    print('upper bound:', float(result['up']))


def check_ye():
    estimator = Estimator(ratio_l2)
    result = estimator.estimate_dilation_regular(
        get_ye_curve(),
        use_face_moments=True, face_dim=0,
        max_depth=5,
        rel_tol_inv=100000,
    )
    result.update(estimator.get_info())
    pprint.pprint(result)
    print('Lower bound is sharp due to the stabilization theorem (Theorem A)')
    print('lower bound:', result['lo'])


def main():
    logging.basicConfig(level=logging.INFO)
    print('The paper by Yuri Malykhin, Evgeny Schepin, "Search of fractal space-filling curves ..."')
    print('The computer-assisted proof of Theorem 1')
    print('')
    print('We consider plain Peano monofractal curves')
    print('')
    print('Genus 2x2:')
    run_search(dim=2, div=2, mult=1, ratio_func=ratio_l2, rel_tol_inv=1000)

    print('Genus 3x3:')
    run_search(dim=2, div=3, mult=1, ratio_func=ratio_l2, rel_tol_inv=1000)

    print('Genus 4x4:')
    run_search(dim=2, div=4, mult=1, ratio_func=ratio_l2, rel_tol_inv=1000)

    print('Genus 6x6')
    print('Plain monofractal curves fall into 3 categories:')
    print('"side", (0,0)->(1,0)')
    print('"diag", (0,0)->(1,1) -- no such curves for even genus')
    print('"med",  (0,0)->(1,1/2)')
    print('Check side curves without diagonal steps:')
    run_search(dim=2, div=6, mult=1, ratio_func=ratio_l2, gate_list=[(SIDE_LINK,)], max_cdist=1, rel_tol_inv=200)
    print('Check median curves without diagonal steps:')
    run_search(dim=2, div=6, mult=1, ratio_func=ratio_l2, gate_list=[(MEDIAN_LINK,)], max_cdist=1, rel_tol_inv=200)

    print('Genus 5x5')
    print('Check diagonal curves without diagonal steps:')
    run_search(dim=2, div=5, mult=1, ratio_func=ratio_l2, gate_list=[(DIAG_LINK,)], max_cdist=1, rel_tol_inv=200)
    print('Check median curves without diagonal steps:')
    run_search(dim=2, div=5, mult=1, ratio_func=ratio_l2, gate_list=[(MEDIAN_LINK,)], max_cdist=1, rel_tol_inv=200)

    print('Check side curves with Non-YE prototype:')
    check_genus5_non_ye()

    print('Check side curves with YE prototype')
    check_genus5_ye_proto()

    print('Check YE curve:')
    check_ye()


if __name__ == "__main__":
    main()
