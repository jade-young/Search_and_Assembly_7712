"""Description: Class for BloomFilterSearchTree object.
"""
import numpy as np
import matplotlib.pyplot as plt
import math

"""Description: Class for BloomFilterSearchTree object which is an collection of BloomFilter objects.
Having a bloom filters arranged into a tree cuts down in the search time as it will return False right
away if a sequence has not been added. It continues traversing the tree in search of the node containing
the sequence if the sequence has been added (with some level of false positivity).
"""
class BloomFilterSearchTree:
    def __init__(self):
        self.bloom_filters = [] # list of bloom filters in the tree
        self.sequences = [] # sequences in the tree

"""Function to add a bloom filter to the tree.
Input:  """