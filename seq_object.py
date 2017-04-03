import seq_generator as sg
import seq_markov as sm
import numpy as np
import bisect
from collections import Counter

uniform_nt_probs = [0.25, 0.25, 0.25, 0.25]
int_to_nt = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
nt_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class sequence():
    def __init__(self, inseq):
        self.dna = inseq
        self.rna = ''
        self.protein = ''
        self.length = len(inseq)

    # def transcribe:

    # def translate:

    # def make_mutation:

    
def generate_seq(nlen, training='ACGT'):
    """Generates a sequence of length nlen using a 0th order Markov model
       with nucleotide probability probs. Returns a string. """

    MC = sm.markov_seq(0, training)
    MC.gen_next_state(nlen)
    return MC.get_chain()
    
x = generate_seq(100)
print x
print Counter(x) # proper proportions?

print sg.filter_zero_prob([1.0/3, 0, 1.0/3, 1.0/3]) # {0: 1.0/3, 2: 1.0/3, 3: 1.0/3}

print sg.transition(sg.filter_zero_prob([1.0/3, 0, 1.0/3, 1.0/3]))

jc = [[0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25]]


def filter_zero_prob(coll):
    idx = xrange(len(coll))
    foo = zip(idx, coll)
    new_coll = [[k, v] for k, v in foo if v != 0]
    return dict(new_coll)


def list_to_dict(coll):
    arr = np.array(coll)
    dimN = arr.shape[0]
    new_dict = {}

    for i in xrange(dimN):
        new_dict.update({i: filter_zero_prob(coll[i])}) #dict(zip(range(dimM), coll[i]))})

    return new_dict
        
    
