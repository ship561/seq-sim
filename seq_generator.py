import random as rand
import numpy as np
from collections import Counter
import bisect
import seq_markov as sm
import seq_object as sobj

print(np.random.poisson(.2,10))

int_to_nt = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
nt_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
uniform_nt_probs = [0.25, 0.25, 0.25, 0.25]

# Jukes-Cantor mutation matrix
jc = [[0, 1.0/3, 1.0/3, 1.0/3],
      [1.0/3, 0, 1.0/3, 1.0/3],
      [1.0/3, 1.0/3, 0, 1.0/3],
      [1.0/3, 1.0/3, 1.0/3, 0]]


def filter_zero_prob(coll):
    idx = xrange(len(coll))
    foo = zip(idx,coll)
    new_coll = [[k, v] for k,v in foo if v != 0]
    return dict(new_coll)


def transition(tr_arr):
    """Determine mutation based on transition probabilities"""
    
    r = rand.random()
    states = tr_arr.keys()
    cum_probs = np.cumsum(tr_arr.values())
    return states[bisect.bisect(cum_probs, r)]
    

# Takes a sequence and mutates it on at an average rate of lambda = mu
def mutate_seq (mu, mut_model, s):
    """Takes a mean number of mutations per sequence, mutation model, and
       a sequence, returns a new sequence that contains mutations
    """
    
    mut_num = np.random.poisson(mu, 1)
    mut_pos = set()
    n = len(s)

    # get the number of mutations to make and select the n positions
    # to mutate
    pos_counter = 0
    while pos_counter < mut_num[0]:
        mut_pos.add(rand.randint(0, n-1))
        pos_counter = len(mut_pos)

    seq = list(s)
    #print mut_num, mut_pos, n, seq
    
    # applies the mutations to the sequence at each mutation position
    for i in mut_pos:
        tr_row = nt_to_int.get(seq[i])
        tr_arr = cum_sum(mut_model[tr_row])
        new_nt = int_to_nt.get(transition(tr_arr))
        seq[i] = new_nt

    return ''.join(seq)


def simulate_n_generations(ngens, mut_model, mu, s):
    """Takes the number of generations (ngen), a mutation model, mean
       number of mutations per sequence, and a sequence and returns a
       boolean for if the sequence ORF was maintained for ngens and
       the final sequence
    """
    
    i = 0
    status = True
    while i < ngens and status:
        s = mutate_seq(mu, mut_model, s)
        i += 1
        if s[:3] != 'ATG' or s[-3:] not in set(['TGA', 'TAG', 'TAA']):
            status = False
            
    # print "completed ", i, "generations", s[:3], s[-3:]
    return [status, s]


def main(nsamp, ngens, mutation_model, mean_mutations, inseq):
    status_counter = {True: 0, False: 0}
    seq_arr = []
    for i in range(nsamp):
        stat, seq = simulate_n_generations(ngens, mutation_model, mean_mutations, inseq)
        status_counter.update({stat: status_counter.get(stat) + 1})
        #seq_arr.append(seq)
    return Counter(status_counter)

#my_seq = 'ATG' + sobj.generate_seq(1000) + 'TGA'

#print mutate_seq(1, jc, my_seq)

#print main(1000, 10000, jc, 10**-9*len(my_seq), my_seq)

def start_state(n, coll):
    if n == 0:
        n = 1
    i = np.random.randint(len(coll)-n-1)
    return coll[i:i+n-1]
    
def build_tr_matrix(n, coll):
    tr_dict = {}

    if n == 0:
        cnts = Counter(coll)
        for k in cnts.keys():
            tr_dict.update({k: cnts})
    else:
        for i in xrange(0, len(coll)-n-1):
            state = tuple(coll[i:i+n+1])
            next_state = coll[i+n+1]
            if state in tr_dict:
                cur = tr_dict.get(state).get(next_state,0) # nested dict call
                tr_dict[state].update({next_state: cur+1})
            else:
                tr_dict[state] = {next_state: 1}

    ## normalize the transition matrix to between [0, 1]
    for state,tr in tr_dict.iteritems():
        next_states = tr.keys()
        cnts = tr.values()
        sigma = sum(cnts)
        props = map(lambda c: float(c)/sigma, cnts)
        tr_dict.update({state: dict(zip(next_states, props))})
            
    return tr_dict

build_tr_matrix(1, list("CATACATGACAGGCTGCTTGGCGAATTCTACGTCAGTACACACCAAGGCTCTGCGCCCGCTGTCGAAAGCGCCTATCGCTAATGTCTGCTGTGGCGCATT"))

def test_tr_n0():
    coll = "CATACATGACAGGCTGCTTGGCGAATTCTACGTCAGTACACACCAAGGCTCTGCGCCCGCTGTCGAAAGCGCCTATCGCTAATGTCTGCTGTGGCGCATT"
    tr_dict = {}
    n = 0
    
    if n == 0:
        cnts = Counter(coll)
        for k in cnts.keys():
            tr_dict.update({k: cnts})

    print tr_dict
