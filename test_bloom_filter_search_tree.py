""" Description: Unit tests for the bloom filter search tree functions and class.
"""
import unittest
from bloom_filter_search_tree import *

class TestInput(unittest.TestCase):
    def test_get_kmers(self):
        expected = ["ACT", "CTG", "TGC", "GCT", "CTA"]
        temp = BloomFilterSearchTree(3, 20, 2)
        observed = temp._get_kmers("ACTGCTA")
        self.assertEqual(expected, observed)
        self.assertEqual(expected, observed)

    def test_check(self):
        temp = BloomFilterSearchTree(3, 20, 2)
        temp.add("ACTGC")
        self.assertTrue(temp.check("CTG"))
        self.assertFalse(temp.check("JADE"))

if __name__ == '__main__':
    unittest.main()