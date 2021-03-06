#!/usr/bin/env python3

import logging
import argparse

from sympy import Rational

import peano.utils as utils
from peano.paths import PathsGenerator
from peano.curves import PathFuzzyCurve
from peano.dilation import Estimator
from peano.gate_utils import GatesGenerator
from peano.subsets import Link


def run_estimator(
        dim, div, pcount,
        finish_max_count=None,
        ratio_func=None, rel_tol_inv=None, rel_tol_inv_mult=None,
        gate_list=None, hyper=False, max_cdist=None,
        upper_bound=None,
        group_by_gates=False,
        output_gates=False, output_stats=False, output_examples=None,
    ):

    if gate_list is not None:
        gates_generator = gate_list
    else:
        gates_generator = GatesGenerator(dim, div, pcount, hyper=hyper).gen_gates()

    def gen_pcurves(gates_iterable):
        for gates_idx, gates in enumerate(gates_iterable):
            if output_gates:
                yield gates
                continue

            logging.info('processing gates: %d', gates_idx + 1)
            paths_gen = PathsGenerator(dim=dim, div=div, links=gates, max_cdist=max_cdist)
            if output_stats:
                counts = []
                for gate in gates:
                    gcnt = len(list(paths_gen.generate_paths_generic(link=gate, std=True)))
                    counts.append(gcnt)
                yield counts
                continue

            kw = {}
            if finish_max_count is not None:
                kw['finish_max_count'] = finish_max_count
            paths_list = list(paths_gen.generate_paths(std=True, **kw))
            logging.warning('gates: %s', [str(g) for g in gates])
            logging.warning('paths: %d', len(paths_list))

            for paths_idx, paths in enumerate(paths_list):
                logging.info('processing gate_paths: %d of %d', paths_idx + 1, len(paths_list))
                yield PathFuzzyCurve.init_from_paths(paths)

    if output_stats:
        scnt = 0
        for counts in gen_pcurves(gates_generator):
            prod = 1
            for cnt in counts:
                prod *= cnt
            print(prod)
            scnt += prod
        print('TOTAL:', scnt)
        return

    if output_gates:
        for gates in gates_generator:
            print(' | '.join([str(gate) for gate in gates]))
        return

    estimator = Estimator(ratio_func, cache_max_size=2**16)

    if group_by_gates:
        pcurve_gens = [(gates, gen_pcurves([gates])) for gates in gates_generator]
    else:
        pcurve_generator = gen_pcurves(gates_generator)
        pcurve_gens = [('all_gates', pcurve_generator)]

    estimate_kwargs = {}
    if rel_tol_inv is not None:
        estimate_kwargs['rel_tol_inv'] = rel_tol_inv
    if rel_tol_inv_mult is not None:
        estimate_kwargs['rel_tol_inv_mult'] = rel_tol_inv_mult

    for gen_id, pcurves_generator in pcurve_gens:
        result = estimator.estimate_dilation_sequence(
            pcurves_generator,
            sat_strategy={'type': 'geometric', 'multiplier': 1.3},
            upper_bound=upper_bound,
            **estimate_kwargs
        )
        print('======')
        print('GENERATOR:', gen_id)
        if not result:
            print('NOT FOUND!')
        else:
            print(result)
            print('lower bound:', float(result['lo']))
            print('upper bound:', float(result['up']))
            if output_examples:
                for curve in result['examples']:
                    print(curve)
                    print('')

        print('', flush=True)


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # configuration args
    argparser.add_argument('--dim', type=int, required=True)
    argparser.add_argument('--pcount', type=int, required=True)
    argparser.add_argument('--div', type=int, required=True)
    argparser.add_argument('--hyper', action='store_true', help='only hypercurves (exit/entrance on hyperfaces)')
    argparser.add_argument('--gates', type=str, help='one tuple of "|"-separated gates')
    argparser.add_argument('--gates-file', type=str)
    argparser.add_argument('--max-cdist', type=int)
    argparser.add_argument('--finish-max-count', type=int)

    # ratio estimation args
    argparser.add_argument('--metric', type=str, choices=['l1','l2','l2_squared','linf'])
    argparser.add_argument('--upper-bound', type=str, help='fractional upper bound, e.g., "11/2"')
    argparser.add_argument('--rel-tol-inv', type=int, help='inverted relative tolerance')
    argparser.add_argument('--rel-tol-inv_mult', type=int, help='multiplier for rel_tol_inv in each epoch')

    # other
    argparser.add_argument('--group-by-gates', action='store_true', help='estimate ratio for each gate')
    argparser.add_argument('--output-gates', action='store_true', help='only gates, do not estimate ratio')
    argparser.add_argument('--output-stats', action='store_true', help='only count, do not work')
    argparser.add_argument('--output-examples', action='store_true', help='print curve examples')
    argparser.add_argument('--verbose', '-v', action='count', default=0, help='loglevel (0=warning, 1=info, 2=debug)')

    args = argparser.parse_args()
    funcs = {
        'l1': utils.ratio_l1,
        'l2': utils.ratio_l2,
        'l2_squared': utils.ratio_l2_squared,
        'linf': utils.ratio_linf,
    }
    if args.dim % 2 == 1 and args.metric == 'l2':
        raise ValueError("Use l2_square for odd dimension!")

    logging.info('args: %s', args)
    gate_list = []
    kwargs = vars(args).copy()
    kwargs['ratio_func'] = funcs.get(kwargs.pop('metric'))
    verbosity = kwargs.pop('verbose')
    if verbosity == 1:
        logging.basicConfig(level=logging.INFO)
    elif verbosity == 2:
        logging.basicConfig(level=logging.DEBUG)

    if (not args.output_stats) and (not args.output_gates) and (args.metric is None):
        raise ValueError("Define metric to estimate ratio!")

    gates_file = kwargs.pop('gates_file')
    if gates_file is not None:
        with open(gates_file) as fh:
            for line in fh:
                gates = [Link.parse_gates(token) for token in line.strip().split('|')]
                gate_list.append(gates)
        kwargs['gate_list'] = gate_list

    gates_str = kwargs.pop('gates')
    if gates_str is not None:
        gates = [Link.parse_gates(token) for token in gates_str.strip().split('|')]
        kwargs['gate_list'] = [gates]

    if args.upper_bound is not None:
        kwargs['upper_bound'] = Rational(args.upper_bound)

    run_estimator(**kwargs)
