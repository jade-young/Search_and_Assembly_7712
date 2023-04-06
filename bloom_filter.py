"""Description: Functions for exploring the optimal array size for the bloom filter as well as
number of hash functions. Class for BloomFilter object.
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
    Input: INT array_size, the size of the array to use as filter; INT num_hash_functions, number of hash functions to use
    Output: BloomFilter object
    '''
    def __init__(self, array_size, num_hash_functions):
        # intialize the array size and the number of hash functions for the new bl0om filter object
        self.array_size = array_size
        self.num_hash_functions = num_hash_functions
        # and set all of the values in the array to False
        self.array = np.zeros(self.array_size, dtype=bool)
        # set the salt values, or seeds (value to be added to the sequence to create a pseudo family of hash functions)
        # bloom filters require that each sequence being added is hashed multiple times to prevent collision
        # the salt value corresponds to each hash function in the "family" used to hash the sequence
        self.salt_values = range(num_hash_functions)

    '''Hash function to get array indices in the bloom filter for this specific sequence
    Input: STRING seq, the sequence being hashed; INT seed, the hash function number
    Output: INTEGER value, the array index
    '''
    def _hash(self, seq, seed):
        return hash(seq + seed) % self.array_size # use modulo so the value is not larger than the array
    
    '''Function to add a sequence to the bloom filter by changing specific indices in the array from 0 to 1
    Input: STRING seq, the sequence being added
    '''
    def add(self, seq):
        for hash_func in self.salt_values:
            hash_result = self._hash(seq, hash_func) # hash the sequence being added
            self.array[hash_result] = 1 # set the element at the hash_result - th index to 1 (the sequence has been added/"seen")
            
    '''Function to check if a sequence has been added to/seen by the bloom filter
    Input: STRING seq, the sequence of interest
    '''
    def check(self, seq):
        for hash_func in self.salt_values:
            hash_result = self._hash(seq, hash_func) # hash the sequence of interest
            if self.array[hash_result] == 0: # check if the element at the hash_result - th index of the array is 0 
                return False # if so, then we can stop and return False because the sequence definitely is not in the bloom filter
        return True # if all of the elements are True then the sequence is (likely) in the bloom filter
