import unittest
import logging
import sys

from quicktions import Fraction

from peano.paths import PathsGenerator
from peano.subsets import Point, Link
from peano.base_maps import BaseMap


class TestGen(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.INFO, stream=sys.stdout)
        self.kws = {'start_max_count': 1, 'finish_max_count': 10 ** 6}

    def test_gen(self):
        configs = [
            {'dim': 2, 'div': 3, 'gates': '(0,0)->(0,1)', 'counts': {'all': 10, 'std': 5, 'hyper': 2, 'hyper_std': 1}},
            {'dim': 2, 'div': 3, 'gates': '(0,0)->(1,1)', 'counts': {'all': 6, 'std': 2, 'hyper': 2, 'hyper_std': 1}},  # hyper:peano
            {'dim': 2, 'div': 4, 'gates': '(0,0)->(0,1)', 'counts': {'all': 298, 'std': 162}},
            {'dim': 2, 'div': 5, 'gates': '(0,0)->(1,1)', 'counts': {'all': 2592, 'std': 659}},
            # 5x5 edge gives 49700 paths, 24850 std; we do not include them here as it slows test

            {'dim': 2, 'div': 3, 'gates': '(0,1/2)->(1,1/2)', 'counts': {'all': 0}},
            {'dim': 2, 'div': 2, 'gates': '(0,0)->(0,1)|(0,0)->(1,1)', 'counts': {'all': 12, 'std': 5}},  # A:4+1*2, B:1*2

            {'dim': 3, 'div': 2, 'gates': '(0,0,0)->(1/2,1/2,1)', 'counts': {'all': 8, 'std': 4}},
        ]
        for conf in configs:
            links = [Link.parse_gates(token) for token in conf['gates'].split('|')]
            gen = PathsGenerator(dim=conf['dim'], div=conf['div'], links=links)
            paths_list = list(gen.generate_paths(**self.kws))
            paths = set(paths_list)
            self.assertEqual(len(paths_list), len(paths))
            self.assertTrue(all(all(path.is_continuous() for path in path_tuple) for path_tuple in paths))

            counts = conf['counts']
            self.assertEqual(len(paths), counts['all'])
            std_paths = set(gen.generate_paths(std=True, **self.kws))
            if 'std' in counts:
                self.assertEqual(len(std_paths), counts['std'])

            if len(links) == 1:
                link = links[0]
                # check that: 1) std are unique 2) give only paths 3) give all paths
                keep_bms = [bm for bm in BaseMap.gen_base_maps(conf['dim']) if bm * link == link and bm != BaseMap.id_map(conf['dim'])]
                seen = set(std_paths)
                for path_tuple in std_paths:
                    path = path_tuple[0]
                    for bm in keep_bms:
                        bm_path = bm * path
                        bm_path_tuple = (bm_path,)
                        self.assertTrue(bm_path == path or (bm_path_tuple not in std_paths))
                        self.assertTrue(bm_path_tuple in paths)
                        seen.add(bm_path_tuple)
                self.assertEqual(len(seen), len(paths))

            if 'hyper' in counts:
                hyper_gen = PathsGenerator(dim=conf['dim'], div=conf['div'], links=links, max_cdist=1)
                hyper_paths = list(hyper_gen.generate_paths(**self.kws))
                self.assertEqual(len(hyper_paths), counts['hyper'])
                std_hyper_paths = list(hyper_gen.generate_paths(std=True, **self.kws))
                self.assertEqual(len(std_hyper_paths), counts['hyper_std'])


    def test_median(self):
        """Get example of median curve"""
        gen = PathsGenerator(dim=2, div=5, links=[Link.parse_gates('(0,0)->(1,1/2)')])
        diag = gen.get_paths_example()
        self.assertIsNotNone(diag)
