""" Description: Unit tests for the bloom filter functions and class.
"""
import unittest
from bloom_filter import *

class TestInput(unittest.TestCase):
    ''' Verify that the calculations are as expected.
    '''
    def test_get_m(self):
        self.assertEqual(get_m(100, 0.01), 959)

    def test_get_k(self):
        self.assertEqual(get_k(100, 959), 7)

    '''Test to make sure the bloom filter is able accurately distinguish from items it has seen 
    and items it hasn't seen... for the most part...
    '''
    def test_bloom_filter_class(self):
        bf = BloomFilter(959, 7)
        bf.add("SEQUENCE")
        bf.add("7712")
        bf.add("DAY3")

        # see if it can recognize the sequences it has seen
        self.assertTrue(bf.check("SEQUENCE"))
        self.assertTrue(bf.check("7712"))
        #self.assertTrue(bf.check("3YAD")) 
        self.assertTrue(bf.check("DAY3"))

        #self.assertFalse(bf.check("7711"))

if __name__ == '__main__':
    unittest.main()