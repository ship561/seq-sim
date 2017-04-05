import re
import seq_markov as mc
from collections import Counter


UNIFORM_NT_PROBS = [0.25, 0.25, 0.25, 0.25]
INT_TO_NT = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
NT_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

JC = {"A": {"C": 1.0/3, "G": 1.0/3, "T": 1.0/3},
      "C": {"A": 1.0/3, "G": 1.0/3, "T": 1.0/3},
      "G": {"A": 1.0/3, "C": 1.0/3, "T": 1.0/3},
      "T": {"A": 1.0/3, "C": 1.0/3, "G": 1.0/3}}


class Sequence():
    """The Sequence object represents elements of a genomic sequence. It
       knows the DNA sequence and thus the corresponding RNA and
       protein sequences as well. Additionally, the sequence object
       can track the mutation rate and a transition matrix to make
       mutations. By default the mutation rate is 10**-9
       errors/site/generation and the transition matrix uses the
       Jukes-Cantor mutation model (i.e. mutation uniform
       probabilities)."""
    
    def __init__(self, inseq, mutation_rate=10**-9, transition_matrix=JC):
        self.dna = inseq
        self.rna = ''
        self.protein = ''
        self.next_seq = ''
        self.length = len(inseq)
        self.functional = True
        self.mu = mutation_rate*self.length
        self.tr = transition_matrix
        
    # def transcribe:

    # def translate:

    def make_mutation(self):
        if self.functional:
            MC = mc.markov_evolution(self.dna, self.tr)
            MC.mutate_seq(self.mu)
            return MC.as_str()
        else:
            return self.dna

    def identify_ORF(self, inseq):
        """Finds the start and stop position of the longest open reading frame
           (ORF) for the DNA sequence
        """
        
        orf_regex = r'(ATG(?:...)*?(?=TAG|TGA|TAA))'
        match_pos = re.finditer(orf_regex, inseq)
        start_pos = 0
        stop_pos = 0
        max_delta = 0
        for m in match_pos:
            tmp_start_pos = m.start()
            tmp_stop_pos = m.end()
            tmp_delta = tmp_stop_pos - tmp_start_pos
            if tmp_delta > max_delta:
                start_pos = tmp_start_pos
                stop_pos = tmp_stop_pos
                max_delta = tmp_delta
        return start_pos, stop_pos, max_delta

    def is_functional(self, orig_seq, mut_seq):
        orig_start, orig_stop, orig_len = self.identify_ORF(orig_seq)
        mut_start, mut_stop, mut_len = self.identify_ORF(mut_seq)
        return orig_start == mut_start \
            and orig_stop == mut_stop \
            and mut_len > 0

    def gen_next_seq(self):
        mut_seq = self.make_mutation()
        if self.is_functional(self.dna, mut_seq):
            self.dna = mut_seq
        else:
            self.functional = False

    def __str__(self):
        return self.dna

    


#x = Sequence('ATG' + generate_seq(99) + 'TGA')
#print x
#print Counter(x.dna) # proper proportions?

#print sg.filter_zero_prob([1.0/3, 0, 1.0/3, 1.0/3]) # {0: 1.0/3, 2: 1.0/3, 3: 1.0/3}

#print sg.transition(sg.filter_zero_prob([1.0/3, 0, 1.0/3, 1.0/3]))
