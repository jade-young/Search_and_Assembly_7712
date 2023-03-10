import numpy as np
import matplotlib.pyplot as plt
import math
# Search...

# EXPLORE
# 1. evaluate the appropriate array size and number of hash functions
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
def evaluate_k(n):
    k_range = ((np.arange(1e6, 1e7, 5e5))/n) * np.log(2)
    plt.plot((np.arange(1e6, 1e7, 5e5)), k_range)
    plt.title(label = "Number of Hash Functions v. Array size")
    #plt.ylim([1e-10, 1e-3])
    #plt.xlim(1, 1e6)
    plt.ylabel("Number of Hash Functions")
    plt.xlabel("Array size")
    plt.savefig('num_hash_array_size_plot.pdf')

def get_m(n, p):
    m = math.ceil((-(n * np.log(p))) / np.log(2)**2)
    return m

def get_k(n, m):
    k = math.ceil((m/n) * np.log(2))
    return k

class BloomFilter:
    def __init__(self, array_size, num_hash_functions):
        self.array_size = array_size
        self.num_hash_functions = num_hash_functions
        self.array = np.zeros(array_size, dtype=np.int8)
    
    def _hash(self, seq, seed):
        value = 0
        for char in seq:
            value = (value * seed) + ord(char)
        return value % self.array_size
    
    def add(self, seq):
        for h in range(self.num_hash_functions):
            result = self._hash(seq, h)
            self.array[result] = 1

    def check(self, seq):
        for h in range(self.num_hash_functions):
            result = self._hash(seq, h)
            if self.array[result] == 0:
                return False
        return True
