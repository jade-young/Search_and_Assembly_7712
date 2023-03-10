import matplotlib.pyplot as plt 
import numpy as np
import time
import random
import math
import reading_input
import pandas as pd
import unittest
from bloom_filter import *

def main():
    time_start = time.time()
    # READ INPUT
    #reads_file_name = "test_fasta.txt"
    reads_file_name = "READS_example.fasta"
    sequence_dict = reading_input.read_fasta(reads_file_name)
    query_file_name = "QUERY_example.fasta"
    query_dict = reading_input.read_query(query_file_name)
    # EXPLORE
    # 1. look at the read count distribution (before error correction)
    reading_input.assess_read_count_distribution(sequence_dict)

    # 2. look at the read length distribution
    reading_input.assess_read_length_distribution(sequence_dict)

    ## CREATE BLOOM FILTER
    # EXPLORE
    # look at different false discovery rates at range of array sizes and number of hash functions
    number_of_items = len(sequence_dict)
    evaluate_m(number_of_items)
    evaluate_k(number_of_items)

    ## get number of sequences
    num_seq = len(sequence_dict)
    ## set desired false positive rate
    false_positive_rate = 0.0001
    ## get optimal array size
    array_size = get_m(num_seq, false_positive_rate)
    ## get optimal number of hash functions
    num_functions = get_k(num_seq, num_seq)

    ## initialize bloom filter object
    bloom_filter = BloomFilter(array_size, num_functions)
    ## add the read sequences to the bloom filter
    for sequence in sequence_dict.values():
        bloom_filter.add(sequence)

if __name__ == '__main__':
    main()