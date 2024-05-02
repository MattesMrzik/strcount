import unittest

import sys
sys.path.append('..')
from genome_str_graph_generator import get_genome_str_graph

class Test_get_genome_str_graph(unittest.TestCase):
    def test_get_genome_str_graph(self):
        segments, links = get_genome_str_graph(
            "resources/config.tsv", "resources/ref.fa", "+", "+", "+", False, True, False)
        true_segments = [["S", "ref_before_prefix_1", "ACGTAACCGGT"],
                         ["S", "ref_prefix_1", "TAAACCCGGGTTT"],
                         ["S", "repeat_1", "AA"],
                         ["S", "ref_suffix_1", "CCCCGGGGT"],
                         ["S", "ref_after_suffix_1", "TTT"]]
        self.assertEqual(segments["ref_1"], true_segments)

        true_links = [["L", "ref_before_prefix_1", "+", "ref_prefix_1", "+", "0M"],
                      ["L", "ref_prefix_1", "+", "repeat_1", "+", "0M"],
                      ["L", "repeat_1", "+", "repeat_1", "+", "0M"],
                      ["L", "repeat_1", "+", "ref_suffix_1", "+", "0M"],
                      ["L", "ref_suffix_1", "+", "ref_after_suffix_1", "+", "0M"]]
        self.assertEqual(links["ref_1"], true_links)


if __name__ == '__main__':
    unittest.main()
