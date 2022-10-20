import unittest
import logging

from peano.gates import GatesGenerator
from peano.subsets import Link


def check_hyper(gate, dim):
    return gate.entrance.face_dim() == dim - 1 and gate.exit.face_dim() == dim - 1


def std_links_list(links_list):
    links_list = [tuple(sorted(link.std() for link in links)) for links in links_list]
    links_list.sort()
    return links_list


class TestGates(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.INFO)

    def test_gen_gates(self):
        configs = [
            {'dim': 2, 'div': 2, 'mult': 1, 'count': 1, 'result': ['(0,0)->(0,1)']},
            {'dim': 2, 'div': 3, 'mult': 1, 'count': 2, 'result': ['(0,0)->(0,1)', '(0,0)->(1,1)']},
            {'dim': 2, 'div': 4, 'mult': 1, 'count': 2, 'result': ['(0,0)->(0,1)', '(0,0)->(1,1/2)']},
            {'dim': 2, 'div': 5, 'mult': 1, 'count': 3, 'result': ['(0,0)->(0,1)', '(0,0)->(1,1)', '(0,0)->(1,1/2)']},

            {'dim': 2, 'div': 2, 'mult': 2, 'only_facet': True,  'count': 1},  # beta-omega
            {'dim': 2, 'div': 3, 'mult': 2, 'only_facet': True,  'count': 2},
            {'dim': 2, 'div': 2, 'mult': 3, 'only_facet': True,  'count': 6},

            {
                'dim': 3, 'div': 2, 'mult': 1,  'count': 6,
                'result': [  # Haverkort's inventory paper
                    '(0,0,0)->(1,0,0)',  # type A
                    '(0,0,0)->(0,1,1)',  # type B
                    '(0,0,0)->(1,1/2,0)',  # type C
                    '(0,0,0)->(1,1/2,1/2)',  # type D
                    '(1/3,0,0)->(1,1/3,1)',  # type E
                    '(0,1/3,1/3)->(2/3,1/3,0)',  # type F
                ],
            },
            {'dim': 3, 'div': 2, 'mult': 1, 'only_facet': True,  'count': 1},
            {'dim': 3, 'div': 2, 'mult': 2, 'only_facet': True,  'count': 35},
        ]
        for conf in configs:
            generator = GatesGenerator(dim=conf['dim'], div=conf['div'], mult=conf['mult'], only_facet=conf.get('only_facet'))
            gen_list = list(generator.gen_gates())
            self.assertEqual(len(gen_list), conf['count'])
            if 'result' in conf:
                true_list = [tuple(Link.parse_gates(token) for token in gates_str.split('|')) for gates_str in conf['result']]
                self.assertEqual(std_links_list(true_list), std_links_list(gen_list))

    def test_hyper_by_all(self):
        for dim, div, mult in [(3, 2, 1)]:  # TODO: add 2,2,2
            all_gates = list(GatesGenerator(dim=dim, div=div, mult=mult).gen_gates())
            hyper_in_all = [gates for gates in all_gates if all(check_hyper(gate, dim) for gate in gates)]
            hyper_in_all.sort()
            hyper_gates = list(GatesGenerator(dim=dim, div=div, mult=mult, only_facet=True).gen_gates())
            hyper_gates.sort()
            self.assertEqual(hyper_in_all, hyper_gates)
