"""Description: Functions for exploring the optimal array size for the bloom filter as well as
number of hash functions. Class for BloomFilter object
"""
import numpy as np
import matplotlib.pyplot as plt
import math
# Portion of search module...

# EXPLORE
# 1. evaluate the appropriate array size and number of hash functions
''' Plot a range of solutions for the equation to find the optimal array size, m
Input: INTEGER n, number of items/sequences that are to be added to the bloom filter
Output: plot showing relationship between false positive rate and array size
'''
def evaluate_m(n):
    # get a range of p (false positive rate)
    p_range = np.arange(1e-10, 1e-3, 1e-10)
    # use the range of p to calculate possible solutions to the array size equation
    m_range = (-(n * np.log(p_range))) / np.log(2)**2
    plt.plot(m_range, p_range)
    plt.title(label = "False Positive Rate v. Array size")
    #plt.ylim([1e-10, 1e-3])
    #plt.xlim(1, 1e6)
    plt.ylabel("False Positive Rate")
    plt.xlabel("Array size")
    plt.savefig('FPR_array_size_plot.pdf')

''' Plot a range of solutions for the equation to find the optimal number of hash functions, k
Input: INTEGER n, number of items/sequences that are to be added to the bloom filter
Output: plot showing relationship between array size and number of has functions
'''
def evaluate_k(n):
    k_range = ((np.arange(1e6, 1e7, 5e5))/n) * np.log(2)
    plt.plot((np.arange(1e6, 1e7, 5e5)), k_range)
    plt.title(label = "Number of Hash Functions v. Array size")
    #plt.ylim([1e-10, 1e-3])
    #plt.xlim(1, 1e6)
    plt.ylabel("Number of Hash Functions")
    plt.xlabel("Array size")
    plt.savefig('num_hash_array_size_plot.pdf')

'''Get the optimal array size
Input: INTEGER n, number of items/sequences that are to be added to the bloom filter;
        FLOAT p, false positive rate
Output: INTEGER m, array size
'''
def get_m(n, p):
    m = math.ceil((-(n * np.log(p))) / np.log(2)**2) # round up to nearest whole number
    return m

'''Get the optimal number of hash functions
Input: INTEGER n, number of items/sequences that are to be added to the bloom filter;
        INTEGER m, array size
Output: INTEGER k, number of hash functions
'''
def get_k(n, m):
    k = math.ceil((m/n) * np.log(2)) # round up to nearest whole number
    return k

''' Class for BloomFilter object
'''
class BloomFilter:
    '''BloomFilter object has a size, number of hash functions, and an array
    Output: BloomFilter ibject
    '''
    def __init__(self, array_size, num_hash_functions):
        self.array_size = array_size
        self.num_hash_functions = num_hash_functions
        self.array = np.zeros(array_size, dtype=np.int8)

    '''Hash function to get array indices in the bloom filter for this specific sequence
    Input: STRING seq, the sequence being hashed; INT seed, the hash function number
    Output: INTEGER value, the array index
    '''
    def _hash(self, seq, seed):
        value = 0
        for char in seq:
            value = (value * seed) + ord(char)
        return value % self.array_size
    
    '''Function to add a sequence to the bloom filter by changing specific indices in the array from 0 to 1
    Input: STRING seq, the sequence being added
    '''
    def add(self, seq):
        for h in range(self.num_hash_functions):
            result = self._hash(seq, h)
            self.array[result] = 1
            
    '''Function to check if a sequence has been added to/seen by the bloom filter
    Input: STRING seq, the sequence of interest
    '''
    def check(self, seq):
        for h in range(self.num_hash_functions):
            result = self._hash(seq, h)
            if self.array[result] == 0:
                return False
        return True
