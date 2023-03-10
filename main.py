# Title: Search and Assembly Driver code
# Description: Uses two files: one FASTA file containing the NGS reads, and another
# containing a query sequence
# 
"""
Title: Search (and Assembly) Driver Code
Author: Jade Young
Date: 03.09.2023
Description: The first section of the program involves reading the input. This input
is expected to be one FASTA file containg the NGS reads, and another contatining
a query sequence. They are stored in dictionaries in which the sequence IDs serve as keys. 
The next part of the program involves creating a bloom filter -- a probabilistic data 
structure that can be used to assess if a given query exists in the set of reads. 
"""

# Import statements
import matplotlib.pyplot as plt 
import numpy as np
import time
import random
import math
import pandas as pd
import unittest
import reading_input # import functions from read_input.py
from bloom_filter import * # import functions from bloom_filter.py

# main function driver code
## TODO: include command line arguments so the end user can input their own files as well as specify their desired false positive rate, output directories, and whether or not to produce the exploration plots
def main():
    # READ INPUT
    #reads_file_name = "test_fasta.txt"
    reads_file_name = "READS_example.fasta" 
    # call the read_reads function which returns the sequencing reads in a dictionary data structure
    sequence_dict = reading_input.read_reads(reads_file_name)
    
    query_file_name = "QUERY_example.fasta"
    # call the read_query function which returns the query sequence(s) in a dictionary data structure
    query_dict = reading_input.read_query(query_file_name)

    # EXPLORE
    # 1. look at the read count distribution
    reading_input.assess_read_count_distribution(sequence_dict) # produces bar plot

    # 2. look at the read length distribution
    reading_input.assess_read_length_distribution(sequence_dict) #produces bar plot

    ## CREATE BLOOM FILTER
    # EXPLORE
    # look at different false discovery rates at range of array sizes and number of hash functions
    num_seq = len(sequence_dict) # get the number of sequences that will be added to the bloom filter
    evaluate_m(num_seq)
    evaluate_k(num_seq)

    # set desired false positive rate
    false_positive_rate = 0.0001
    # get optimal array size
    array_size = get_m(num_seq, false_positive_rate)
    # get optimal number of hash functions
    num_functions = get_k(num_seq, num_seq)

    # initialize bloom filter object
    bloom_filter = BloomFilter(array_size, num_functions)
    # add the read sequences to the bloom filter
    for sequence in sequence_dict.values():
        bloom_filter.add(sequence)


if __name__ == '__main__':
    main()