"""Description: Class for BloomFilterSearchTree object.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import bisect
from bloom_filter import * # import functions from bloom_filter.py

"""Description: Class for BloomFilterSearchTree object which is an collection of BloomFilter objects.
Having a bloom filters arranged into a tree cuts down in the search time as it will return False right
away if a sequence has not been added. It continues traversing the tree in search of the node containing
the sequence if the sequence has been added (with some level of false positivity).
"""
class BloomFilterSearchTree:
    """Initialize the bloom filter tree object.
    Input: INT k, the k-mer length; INT array_size, the size of the array to be used as a bloom filter;
            INT num_hash_functions, the number of hash functions to use for the bloom filters
    Output: BloomFilterSearchTree object"""
    def __init__(self, k, array_size, num_hash_functions):
        self.k = k
        self.array_size = array_size
        self.num_hash_functions = num_hash_functions
        # dictionary to hold the bloom filters where each key is a node in the tree representing a sequence
        # and its value is the corresponding bloom filter for all of the k-mers of that sequence
        self.bloom_filter_search_tree = {} 
    
    """Function to generate all of the k-mers for a given sequence. Each of these k-mers will then be added to a bloom filter
    for the sequence.
    Input: STRING seq, the sequence being added to the bloom filter tree
    Ouput: LIST kmers, a list of the possible kmers for the given sequence"""
    def _get_kmers(self, seq):
        kmers = []
        for i in range(len(seq) - self.k +1):
            kmers.append(seq[i:i+self.k])
        return kmers

    """Function to add a bloom filter to the tree. This is done by creating a bloom filter for the input sequence
    to which each of its kmers are added. The given sequence becomes a node and is stored as a key in the 
    BloomFilterSearchTree object's dictionary, and the bloom filter is the corresponding value.
    Input: STRING seq, sequence being added"""
    def add(self, seq):
        # first get the kmers for the sequence
        seq_kmers = self._get_kmers(seq)
        # construct nodes by overlapping the sequence kmers
        # pre/suffix overlap the kmer by one base to aid in the construction of the de bruijn graph
        for i in range(len(seq_kmers)):
            if i == 0 : # if it is the first kmer
                prefix == "" #then don't add anything to the front of the kmer
            else: 
                prefix = seq_kmers[i-1][-self.k+1:]
            if i == len(seq_kmers) - 1: # if it is the last kmer
                suffix = "" #then don't add anything to the end of the kmer
            else:
                suffix = seq_kmers[i+1][:self.k-1]
            node = prefix + seq_kmers[i] + suffix
            if node not in self.bloom_filter_search_tree: 
                self.bloom_filter_search_tree[node] = BloomFilter(self.array_size, self.num_hash_functions)
            self.bloom_filter_search_tree[node].add(seq_kmers[i])

    """Function to check if a sequence is in the bloom filter search tree
    Input: STRING seq, the sequence of interest/query
    Output: BOOL eval, True if the sequence is (likely) in the tree and False if it is not """
    def check(self, seq):
        # loop through the bloom filters
        for pointer, bloom_filter in enumerate(self.bloom_filters):
            if bloom_filter.check(seq): # call the bloom filter's check function
                index = bisect.bisect_left()
