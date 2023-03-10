"""Description: Functions related to reading input as well as exploring the read 
count and length distributions.
"""

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd

# module 1: read in fasta file and query file
'''
Read in FASTA file contatining the sequencing reads
Input: STRING reads_file_name
Output: DICTIONARY sequence_dictionary, where key is a sequence ID and value is corresponding sequence
'''
def read_reads(reads_file_name):
    sequence_dictionary = {}
    with open(reads_file_name) as file:
        while True:
            sequence_ID = file.readline().strip()[1:]
            sequence = file.readline().strip()
            if sequence_ID != "": # don't add an empty string to the dictionary
                sequence_dictionary[sequence_ID] = sequence # break when there aren't anymore seqeunces in the file
            if not sequence: break
    return sequence_dictionary

'''
Read in FASTA file contatining the initial query
Input: STRING query_file_name
Output: DICTIONARY query_dictionary, where key is a sequence ID and value is corresponding sequence
'''
def read_query(query_file_name):
    query_dictionary = {}
    with open(query_file_name) as file:
        while True:
            query_ID = file.readline().strip()[1:]
            query = file.readline().strip()
            if query_ID != "": # don't add an empty string to the dictionary
                query_dictionary[query_ID] = query # break when there aren't anymore sequences in the file
            if not query: break
    return query_dictionary

# EXPLORE
# 1. report the count distribution 
# could also help to choose a threshold for the filtering step during error correction later
'''
Record the frequency of each read and store in a dictionary
Then plot the read count on a barplot
Input: DICTIONARY read_dict (output from read_reads function)
Output: bar plot
'''
def assess_read_count_distribution(read_dict):
    # first get a dictionary that records the number of times each unique sequence appears in the reads dictionary
    read_count = {} # where each key is a read and the value is the count (MAYBE THE KEY SHOULD BE THE ID?)
    for read_ID in read_dict.keys():
        sequence = read_dict[read_ID]
        if sequence in read_count:
            read_count[sequence] +=1
        else:
            read_count[sequence] = 1

    # look at the read count distribution (before error correction if applicable)
    # sort the dictionary entries by ascending read count
    read_count_dict_sorted = dict(sorted(read_count.items(),key=lambda x:x[1]))
    # plot the read count distribution
    ##plt.bar(list(read_count_dict_sorted.keys()), read_count_dict_sorted.values(), 1.0, color='g', tick_label = None)
    read_count_sorted_list = [key for key, val in read_count_dict_sorted.items() for _ in range(val)] # need a list for the hist function
    df = pd.DataFrame(list(read_count_dict_sorted.items()), columns=['Sequence', 'Frequency'])
    plot = df.plot(kind='bar', x='Sequence', ylabel="Frequency", legend = False, title="Sequence Frequency")
    plot.tick_params(axis="x", which="both", top = False, bottom = False, labelbottom=False)
    fig = plot.get_figure()
    fig.savefig('read_count_dist_plot.pdf')
    #plt.hist(read_count_sorted_list, bins=100)
    ##plt.tick_params(axis="x", which="both", top = False, bottom = False, labelbottom=False)
    #plt.ylabel("Read Frequency")
    #plt.axvline(np.mean(list(read_count_dict_sorted.values())), color='k', linestyle='dashed', linewidth=1)
    #plt.savefig('read_count_distribution_plot.pdf')  
    

# 2. report read length distribution 
# could help to also establish a read length threshold to further filter the reads
'''
Record the frequency of each read unique read lengthand store in a dictionary
Then plot the read length on a barplot.
Input: DICTIONARY read_dict (output from read_reads function)
Output: bar plot
'''
def assess_read_length_distribution(read_dict):
    read_length_dict = {}
    for key in read_dict.keys():
        read_length = len(read_dict[key])
        if read_length in read_length_dict:
            read_length_dict[read_length] += 1
        else:
            read_length_dict[read_length] = 1
    
    # look at the read length distribution (before error correction if applicable)
    # sort the dictionary entries by read length
    read_length_dict_sorted = dict(sorted(read_length_dict.items(),key=lambda x:x[1], reverse = True))
    
    #read_length_sorted_list = [val for key, val in read_length_dict_sorted.items() for _ in range(val)]

    #plt.hist(read_length_sorted_list, bins=100)
    #plt.bar(list(read_length_dict_sorted.keys()), read_length_dict_sorted.values(), 1.0, color='g', tick_label = None)
    df = pd.DataFrame(list(read_length_dict_sorted.items()), columns=['Read Length', 'Frequency'])
    plot = df.plot(kind='bar', x='Read Length', ylabel="Frequency", legend = False, title="Read Length Frequency")
    fig = plot.get_figure()
    fig.savefig('read_length_dist_plot.pdf')


    ##plt.tick_params(axis="x", which="both", top = False, bottom = False, labelbottom=False)
    #plt.ylabel("Read Length")
    #plt.savefig('read_length_distribution_plot.pdf') 

