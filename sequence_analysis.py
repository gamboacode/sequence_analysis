import numpy as np
import random
import re
from itertools import product
import math

"""In this module, we provide functions and classes that are used in 
DNA sequence analysis, with the final goal to implement High order 
Markov Chains with context."""


def dna_to_rna(s):
    """Input:
            DNA string s
            
        Output:
            RNA string or ValueError if some character in s is 
            not in "ACGT" """
    rna = ""
    for nucleotide in s:
        if nucleotide not in "ACGT":
            raise ValueError("Invalid character: String s must have only characters A, C G or T")
        elif nucleotide == "T":
            rna += "U"
        else:
            rna += nucleotide
    return rna


def get_dict():
    """Reads the hash table for the (rna sequence, amino acid) pair and outputs 
        the corresponding pair in a dictionary where the rna sequence is the key
        and the value is the amino acid"""
    dic = {}
    with open("table.txt", "r") as my_file:
        for line in my_file:
            pattern_line = re.findall(r'([A-Z]{3,3}) ([A-Z]|\*)', line)
            if len(pattern_line) != 4:  # Each line should have 4 pairs
                continue
            else:
                for pattern_tuple in pattern_line:
                    dic[pattern_tuple[0]] = pattern_tuple[1]
    return dic


def get_dict_list():
    """Reads the hash table for the (codon, amino acid) pair and outputs 
        the corresponding pair in a dictionary where the codon is the key"""
    
    dic = {}
    with open("table.txt", "r") as my_file:
        for line in my_file:
            pattern_line = re.findall(r'([A-Z]{3,3}) ([A-Z]|\*)', line)
            if len(pattern_line) != 4: #Each line should have 4 pairs
                continue
            else:
                for pattern_tuple in pattern_line:
                    if pattern_tuple[1] not in dic:
                        dic[pattern_tuple[1]] = []
                    dic[pattern_tuple[1]].append(pattern_tuple[0])
    return dic


def rna_to_prot(s):
    """Input:
            RNA or DNA string s
            
        Output:
            Aminoacid string or ValueError if some character in s is 
            not in "ACGU" """
    s = dna_to_rna(s)  # Convert DNA to RNA if needed
    dict_codon_to_amino = get_dict()
    amino = ""
    if len(s) % 3 != 0:
        raise ValueError("Invalid length: String s must have a length that is a multiple of 3")
    for i in range(0, len(s), 3):
        codon = s[i:i+3]
        if codon not in dict_codon_to_amino:
            raise ValueError("Invalid codon: String s must have only \\\
                             characters A, C, G or U")
        amino += dict_codon_to_amino[codon]
    return amino


def get_probability_dict():
    """
    Return the probability dictionary with codons as keys and probabilities as values.
    The probabilities are taken over the space of codons that represent an amino acid.
    """
    dic_init = {amino:[] for amino in get_dict_list()}
    with open("table.txt", "r") as my_file:
        for line in my_file:
            patt_line = re.findall(r'([A-Z]{3,3}) ([A-Z]|\*).{11}(\(\s?\d*\))', line)
            if len(patt_line) != 4:
                continue
            else:
                for patt_tuple in patt_line:
                    amino = patt_tuple[1]
                    dic_init[amino].append((patt_tuple[0], int(patt_tuple[2].strip("() "))))
    
    dic_total = {}
    for amino in dic_init:
        total = 0
        for tuple in dic_init[amino]:
            total += tuple[1]
        dic_total[amino] = total
    
    dic_codons = {}
    for amino in dic_init:
        for tuple in dic_init[amino]:
            dic_codons[tuple[0]] = tuple[1] / dic_total[amino]
    
    return dic_codons



class ProteinToMaxRNA:
    """This class converts a protein sequence into the most likely RNA sequence
    that could have produced it, based on the codon probabilities."""
    prob_dic = get_probability_dict()
    aa_to_codon = get_dict_list()
    
    def __init__(self):
        pass

    def convert(self, s):
        """Converts a protein sequence s into the most likely RNA sequence
        to be the source of the protein sequence s."""
        rna_seq = []
        for aa in s:
            codons = ProteinToMaxRNA.aa_to_codon[aa]
            most_likely = codons[0]
            for codon in codons:
                if ProteinToMaxRNA.prob_dic[codon] > ProteinToMaxRNA.prob_dic[most_likely]:
                    most_likely = codon 
            rna_seq.append(most_likely)
        return "".join(rna_seq)


def random_event(dist):
    """
    Takes as input a dictionary from events to their probabilities.
    Return a random event sampled according to the given distribution.
    The probabilities must sum to 1.0
    """
    out = random.uniform(0, 1)
    ordered_dic = list(dist.items())
    dic_prob = {}

    init = 0.0
    for tuple in ordered_dic:
        dic_prob[tuple[0]] = [init, init + tuple[1]]
        init += tuple[1]

    for event in dic_prob:
        if dic_prob[event][0] <=  out <= dic_prob[event][1]:
            return event
        

class ProteinToRandomRNA(object):
    """This class converts a protein sequence into a random RNA sequence
    according to the codon probability distribution."""
    prob_dic = get_probability_dict()
    aa_to_codon = get_dict_list()
    
    def __init__(self):
        pass

    def convert(self, s):
        """Converts the input protein sequence s into a random 
        RNA sequence according to the propability distribution"""
        rna_seq = []
        for aa in s:
            codons = ProteinToRandomRNA.aa_to_codon[aa]
            probabilities = [ProteinToRandomRNA.prob_dic[codon] for codon in codons]
            distribution = dict(zip(codons, probabilities))

            rna_seq.append(random_event(distribution))
        
        return "".join(rna_seq)
    

def sliding_window(s, k):
    """
    Returns a generator that can be iterated over all starting positions
    of a k-window in the sequence. For each starting position the generator
    returns the nucleotide frequencies in the window as a dictionary.

    Example:
    >>> s = "ACGTACGT"
    >>> for freq in sliding_window(s, 3):
    ...     print(freq)
    {'A': 1, 'C': 1, 'T': 0, 'G': 1}
    {'A': 0, 'C': 1, 'T': 1, 'G': 1}
    {'A': 1, 'C': 0, 'T': 1, 'G': 1}
    {'A': 1, 'C': 1, 'T': 1, 'G': 0}
    {'A': 1, 'C': 1, 'T': 0, 'G': 1}
    {'A': 0, 'C': 1, 'T': 1, 'G': 1}
    {'A': 0, 'C': 0, 'T': 2, 'G': 1}
    """
    start_index = 0
    end_index = k
    for _ in s:
        if end_index <= len(s):
            subword = s[start_index:end_index]
        else:
            break
        
        bases = ["A", "C", "T", "G"]
        dic = { base: subword.count(base) for base in bases}
        yield dic

        start_index += 1
        end_index += 1


def context_list(s, k):
    """Returns dictionary where the keys are k-mers and the 
    values are the concatenation of symbols appearing after each k-mer.
    
    Example:
    >>> s = "AACGTAACGT"
    >>> k = 3
    >>> context_list(s, k)
    {'AAC': 'GG', 'ACG': 'T', 'CGT': 'A', 'GTA': 'A', 'TAA': 'C', 'ACG': 'T'}
    """
    start = 0
    dic = {}
    while start + k <= len(s) - 1:
        kmer = s[start:start + k]
        if kmer in dic:
            dic[kmer].append(s[start + k])
        else:
            dic[kmer] = [s[start + k]]
        start += 1
    for key, val in dic.items():
        dic[key] = "".join(val)
    return dic


def context_probabilities(s, k):
    """Converts frequencies from the dictionary context_list(s, k)
    to probabilities for each k-mer
    """
    context_dic = context_list(s, k)
    dic = {}
    proteins = list("ACTG")
    for kmer in context_dic:
        dic[kmer] = {protein: (context_dic[kmer].count(protein) / len(context_dic[kmer])) for protein in proteins}
    return dic


def better_context_probabilities(s, k):
    start = 0
    dic = {}
    proteins = list("ACTG")
    while start + k < len(s) - 1:
        kmer = s[start:start + k]
        if kmer not in dic:
            dic[kmer] = {protein:0.0 for protein in proteins}
        dic[kmer][s[start + k]] += 1.0
        start += 1
    
    dic_probs = {}
    for kmer, freq in dic.items():
        total = sum(freq.values())
        dic_probs[kmer] = {protein: freq[protein] / total for protein in proteins}

    for kmer in product("ACGT", repeat=k):
        kmer_joint = "".join(kmer)
        if kmer_joint not in dic_probs:
            dic_probs[kmer_joint] = dict(zip("ACGT", [0.25, 0.25 , 0.25, 0.25]))
    
    return dic_probs


def context_pseudo_probabilities(s, k):
    """Completes the context probabilities by adding pseudo-counts
    to ensure that all k-mers have a probability distribution."""
    context_prob = better_context_probabilities(s, k)
    all_kmers = ["".join(tuple) for tuple in product("ACGT", repeat=k)]

    for kmer in all_kmers:
        if kmer not in context_prob:
            context_prob[kmer] = dict(zip("ACGT", [0.25 for _ in range(4)]))

    return context_prob



class MarkovChain:    
    """This class implements a Markov Chain of kth order for RNA sequences.
    It uses zeroth order probabilities for the first k characters and the kth order
    probabilities for the subsequent characters."""

    def __init__(self, zeroth, kth, k=2):
        """Initializes the Markov Chain of kth order"""
        self.k = k
        self.zeroth = zeroth
        if len(kth) != 4 ** k:
            raise ValueError(f"Invalid kth order probabilities: expected {4 ** k} k-mers, got {len(kth)}")
        self.kth = kth
        
    def generate(self, n, seed=None):
        """Generates a random RNA sequence of length n using the 
        Markov Chain defined by zeroth and kth order probabilities."""

        if seed is not None:
            np.random.seed(seed=seed)
            random.seed(a=seed)
        
        left = n-self.k
        if left < 0:
            return
        
        rna_seq = ""
        for _ in range(self.k):
            protein = random_event(self.zeroth)
            rna_seq += protein # type: ignore
        
        for _ in range(left):
            kmer = rna_seq[-self.k:]
            protein = random_event(self.kth[kmer])
            rna_seq += protein # type: ignore
        
        return rna_seq

    def probability(self, s):
        """Calculates the probability of a given RNA sequence s
        under the Markov Chain defined by zeroth and kth order probabilities."""
        prob = 1.0

        for i in range(self.k):
            prob *= self.zeroth[s[i]]

        for i in range(self.k, len(s)):
            prob *= self.kth[s[i - self.k:i]][s[i]]

        return prob
    
    def log_probability(self, s):
        """Calculates the log probability of a given RNA sequence s
        under the Markov Chain defined by zeroth and kth order probabilities.w"""
        log_prob = 0.0

        for i in range(self.k):
            log_prob += math.log2(self.zeroth[s[i]])

        for i in range(self.k, len(s)):
            log_prob += math.log2(self.kth[s[i - self.k:i]][s[i]])

        return log_prob
    

class SimpleMarkovChain(object):
    """This class implements a simple Markov Chain of kth order for RNA sequences.
    It uses the k-th context of the input sequence to generate new sequences."""
    def __init__(self, s, k):
        self.k = k
        self.s = s

    def generate(self, n, seed=None):
        """Generates a random RNA sequence of length n using a simple Markov Chain
        based on the k-th context of the input sequence s."""
        np.random.seed(seed=seed)
        
        left = n-self.k

        if left < 0:
            return
        
        context_dic = context_list(self.s, self.k)

        rna_seq = ""
        for i in range(self.k):
            protein = np.random.choice(np.array(list("ACGT")))
            rna_seq += protein
        
        for i in range(left):
            kmer = rna_seq[-self.k:]
            if kmer in context_dic:
                protein = np.random.choice(np.array(list(context_dic[kmer])))
            else:
                protein = np.random.choice(np.array(list("ACGT")))
            rna_seq += protein
        
        return rna_seq


def kmer_index(s, k):
    """Returns a dictionary where the keys are k-mers and the values are lists of
    starting indices of the k-mers in the string s.
    
    Example:
    >>> s = "AACGTAACGT"
    >>> k = 3
    >>> kmer_index(s, k)
    {'AAC': [0, 5], 'ACG': [1, 6], 'CGT': [2, 7], 'GTA': [3], 'TAA': [4]}
    """
    dic = {}
    for i in range(len(s) - k + 1):
        kmer = s[i:i+k]
        if kmer not in dic:
            dic[kmer] = [i]
        else:
            dic[kmer].append(i)
    return dic


def codon_probabilities(rna):
    """
    Given an RNA sequence, simply calculates the probability of
    all 3-mers empirically based on the sequence
    """
    dic = {"".join(codon): 0.0 for codon in product("ACGU", repeat=3)}
    for i in range(len(rna)-2):
        three_mer = rna[i:i+3]
        dic[three_mer] += 1.0

    for three_mer in dic:
        dic[three_mer] *= 1/(len(rna) - 2)

    return dic

def kullback_leibler(p, q):
    """Given two probability distributions p and q, calculates the 
    Kullback-Leibler divergence D(p || q). The distributions p and q 
    should be dictionaries where keys are events and values are their 
    probabilities."""
    suma = 0
    for event in p:
        if q[event] != 0 or p[event] == 0:
            suma += p[event] * math.log2(p[event] / q[event])
    return suma


def get_stationary_distributions(transition):
    """
    The function get a transition matrix of a degree one Markov chain as parameter.
    It returns a list of stationary distributions, in vector form, for that chain.
    The stationary distribution is a vector that satisfies the equation πP = π,
    where P is the transition matrix and π is the stationary distribution.
    """
    transpose = transition.T
    eigenvalues, eigenvectors = np.linalg.eig(transpose)

    stationary_vectors = eigenvectors.T[eigenvalues == 1.0]
    for vector in stationary_vectors:
        vector *= (1/np.sum(vector))
    
    return stationary_vectors