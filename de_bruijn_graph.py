"""de Bruijn Graph class and associated functions. Constructs a de bruijn graph from the reads and finds the longest
contig containing the query sequence.
NOTE: this class is only used when a query sequece is determined to be in the Bloom Filter Search Tree is found."""

from collections import defaultdict # this method is used to define the type for the values in the dictionaries used

class DeBruijnGraph:
    """Initialize the de Bruijn Graph object.
    Input: DICT reads, dictionary of sequencing reads for which the key is the sequence ID and the value is the
    sequence string; INT k, length of the k-mers
    Output:
    """
    def __init__(self, reads, k):
        self.k = k
        # build the k-mer table that will be used to construct the graph
        self.kmer_table = self._build_kmer_table(self.k, reads)
        # construct the graph
        self.graph = self._construct_debruijn_graph(self.kmer_table)

    
    """Build a k-mer table. All of the possible k-mers are configured for each sequencing read
    and its frequency is recorded. 
    Input: INT k, length of the k-mers; DICT reads, dictionary of sequenceing reads where the key is the 
    sequence ID and the value is a sequence string
    Ouptut: DICT kmer_table, dictionary where each key is a kmer and the value is a count for that k-mer"""
    def _build_kmer_table(self, k, reads):
        kmer_table = defaultdict(int)
        for seq in reads.values():
            for i in range(len(seq) - k+ 1):
                kmer = seq[i:i+k]
                kmer_table[kmer] +=1
        return kmer_table

    """Construct the de Bruijn Graph. Returned is a dictionary of k-1 mer prefixes (the keys) and
    their corresponding suffixes (the values). Following the overlap when traversing the tree will help to
    find the longest contig containing the query sequence in the get_longest_contig method below.
    Input: DICT kmer_table, dictionary of k-mers from all of the sequencing reads
    Output: DICT de_bruijn_graph, where the key is the prefix of a kmer and the value is a list of possible suffixes."""
    def _construct_debruijn_graph(self, kmer_table):
        de_bruijn_graph = {}
        for kmer in kmer_table:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            if prefix not in de_bruijn_graph:
                #de_bruijn_graph[prefix][suffix] = count
                de_bruijn_graph[prefix] = [suffix]
            else:
                de_bruijn_graph[prefix].append(suffix)
        return de_bruijn_graph

    """Construct the longest contig that contains the query sequence
    Input:
    Output:"""
    def get_longest_contig(self, query_seq):
        # start with the first k-mer of the query
        first_kmer = query_seq[:self.k]
        #print(first_kmer)
        # then split into left and right
        left_kmer = first_kmer[:-1]
        right_kmer = first_kmer[1:]
        #print(left_kmer)
        #print(right_kmer)
        # start at the left k-mer
        # keep track of the nodes visited and th edges between them
        contig_nodes_left = []
        #print(contig_nodes_left)
        contig_edges = []
        current_node = left_kmer
        while True:
            # get the k-1 mers that overlap with the current node
            neighbors = self.graph[current_node]
            #print(neighbors)
            if len(neighbors) == 1:
                # if there is only one neighbor, set it as the next node and loop
                next_node, edge_label = neighbors[0]
                contig_nodes_left.append(current_node)
                contig_edges.append(edge_label)
                current_node = next_node
            else:
                # branching point or dead end...?
                contig_nodes_left.append(current_node)
                break
            
        # do the same for the right kmer
        contig_nodes_right = []
        current_node = right_kmer
        #print(current_node)
        while True:
            neighbors = self.graph[current_node]
            #print(neighbors)
            if len(neighbors) == 1:
                next_node, edge_label = neighbors[0]
                contig_nodes_right.append(current_node)
                contig_edges.append(edge_label)
                current_node = next_node
            else:
                # branching point/dead end...
                contig_nodes_right.append(current_node)
                break
        
        # reverse the order of the right nodes
        contig_nodes_right.reverse()
        contig = contig_nodes_left[:-1] + first_kmer + contig_nodes_right[1:]

        # join the nodes with the query to get the final contig
        return ''.join(contig)


