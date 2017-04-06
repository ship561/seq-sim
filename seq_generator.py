from collections import Counter
import seq_markov as mc
import seq_object as sobj
import sys


def generate_seq(nlen, training='ACGT'):
    """Generates a sequence of length nlen using a 0th order Markov model
       with nucleotide probability probs. Returns a string. """

    MC = mc.markov_seq(0, training)
    MC.gen_next_state(nlen)
    return MC.get_chain()


def simulate_n_generations(ngens, s):
    """Takes the number of generations (ngen), a mutation model, mean
       number of mutations per sequence, and a sequence and returns a
       boolean for if the sequence ORF was maintained for ngens and
       the final sequence
    """
    
    i = 0
    while i < ngens and s.functional:
        s.make_mutation()
        i += 1
            
    return s


def main(nsamp, ngens, mutation_model, mutation_rate, inseq):
    status_counter = {True: 0, False: 0}
    for _ in xrange(nsamp):
        seq_obj = sobj.Sequence(inseq, mutation_rate, mutation_model)
        seq = simulate_n_generations(ngens, seq_obj)
        status_counter.update({seq.functional: status_counter.get(seq.functional) + 1})
        #seq_arr.append(seq)
    return Counter(status_counter)

x = 'ATG' + generate_seq(99) + 'TGA'
print x
print Counter(x)

print main(10, 1000, sobj.JC, 10**-5, x)


# sys.argv[1]

jc = [[0, 1.0/3, 1.0/3, 1.0/3],
      [1.0/3, 0, 1.0/3, 1.0/3],
      [1.0/3, 1.0/3, 0, 1.0/3],
      [1.0/3, 1.0/3, 1.0/3, 0]]

isinstance(jc, list)
isinstance(np.array(jc), np.ndarray)
isinstance({'a': 1}, dict)

def filter_zero_prob(coll):
    """Special case for filtering out the zero probability transitions
    from the matrix"""
    
    idx = xrange(len(coll))
    foo = zip(idx, coll)
    new_coll = [[k, v] for k, v in foo if v != 0]
    return dict(new_coll)

def list_to_dict(coll):
    """takes a list and makes it a dictionary"""
    
    arr = np.array(coll)
    dimN = arr.shape[0]
    new_dict = {}
    
    for i in xrange(dimN):
        new_dict.update({i: self.filter_zero_prob(coll[i])})
        
    return new_dict

