""" Description: Unit tests for the de Bruijn Graph class.
"""
import unittest
from de_bruijn_graph import *

class TestInput(unittest.TestCase):
    def test_build_kmer_table(self):
        reads = {
            "SEQ1": "ATGCTA",
            "SEQ2": "GCTAGC",
            "SEQ3": "TAGCAC",
            "SEQ4": "GCACAT",
            "SEQ5": "ACATGC"
        }
        k = 3
        temp = DeBruijnGraph(reads, k)
        expected_kmer_table = {
            "ATG": 2,
            "TGC": 2,
            "GCT": 2,
            "CTA": 2,
            "TAG": 2,
            "AGC": 2,
            "GCA": 2,
            "CAC": 2,
            "ACA": 2,
            "CAT": 2
        }
        #expected_graph = {
        #    "AT": {"TG": 2},
        #    "TG": {"GC": 2},
        #    "GC": {"CT": 2},
        #    "CT": {"TA": 2},
        #    "TA": {"AG": 2},
        #    "AG": {"GC": 2},
        #    "GC": {"CA": 2},
        #    "CA": {"AC": 2},
        #    "AC": {"CA": 2},
        #    "CA": {"AT": 2}
        #}
        
        expected_graph = {
            "AT": ["TG"],
            "TG": ["GC"],
            "GC": ["CT", "CA"],
            "CT": ["TA"],
            "TA": ["AG"],
            "AG": ["GC"],
            "GC": ["CA"],
            "CA": ["AC"],
            "AC": ["CA"],
            "CA": ["AC", "AT"]
        }
        observed_kmer_table = temp.kmer_table
        observed_graph = temp.graph
        self.assertEqual(expected_kmer_table, observed_kmer_table)
        self.assertEqual(expected_graph, observed_graph)


if __name__ == '__main__':
    unittest.main()