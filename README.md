# Search And Assembly
## Goal
Given a FASTA file containing a set of sequencing reads and a file containing a query sequence, generate the largest contig that can be constructed from the reads that contain the query sequence. 
## Description
Researchers may want to search a genome for a sequence of interest to verify a knockin/out experiment, relate one organism to another based on the shared presence (or absence) of a specific gene, to identify regulatory regions surrounding a gene to name a few scenarios. However, if they are studying an organism that is not widely studied and a reference genome does not exist, one must assemble the genome using only the sequencing reads to perform the aforementioned searches. This is called de novo genome assembly. During sequencing, the DNA of the organism is fragmented and each short fragment is sequenced. The copies of these fragments, or sequencing reads, need to be arranged in the correct order (and orientation) to get back the original sequence. One way to do this is by constructing a de Bruijn graph in which each node is a sequence of length k. When each node is visited once, the graph is able to reconstruct long stretches of the sequence. A bloom filter can be used to search the reads for the query sequence.

## Installation Instructions
Download all files to local computer. 
### Dependencies
This program is written in Python version 3.10.8, and has the following dependencies:

+ matplotlib.pyplot
+ numpy
+ time
+ random
+ math
+ pandas
+ unittest

If any of the above libraries are not installed, install them using `pip.` For instance: `pip install numpy`

## Examples of Use
To run this program, change directories (`cd`) into the folder containing the code and input files.
Once in the directory, run using the following command: `python3 main.py`
## Input and Output Descriptions
_**NOTE**: All input is expected to be in the same folder as the Python files. Any output generated will be added to the same directory in which `python3 main.py` is run._
### Input
Two input files are expected:

+ a FASTA file containing the sequencing reads
+ a FASTA file containing the query sequence(s)

The expected format for both is standard FASTA format includes a comment line starting with '>' followed by the sequence on the next line.
For example:
```
>sequence_ID_1
YYYYYYYYYYYYYYYYYYYYYYYYYYYY
>sequence_ID_2
ZZZZZZZ
```
### Output
The current output should be 4 figures:

+ a barplot showing the read count frequency
+ a barplot showing the read length frequency
+ a plot showing the relationship between false positive rates and array sizes
+ a plot showing the relationship between array sizes and number of hash functions
