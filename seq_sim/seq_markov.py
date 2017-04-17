import numpy as np
from collections import Counter
import bisect
import seq_utils as su


class markov_chain():
    def __init__(self, order=0, start_state=None, transition_matrix=None):
        """Generic Markov Chain object. This object is meant to be
           generalizable to any sort of Markov process. It only needs
           to know how to transition to the next state.
        """
        self.order = order
        self.start_state = start_state
        self.transition_matrix = transition_matrix

    def make_transition(self, tr_arr):
        """Determine mutation based on transition probabilities"""

        r = np.random.rand()
        states = tr_arr.keys()
        cum_probs = np.cumsum(tr_arr.values())
        return states[bisect.bisect(cum_probs, r)]

    def filter_zero_prob(self, coll):
        """Special case for filtering out the zero probability transitions
           from the matrix"""

        idx = xrange(len(coll))
        foo = zip(idx, coll)
        new_coll = [[k, v] for k, v in foo if v != 0]
        return dict(new_coll)

    def list_to_dict(self, coll):
        """takes a list and makes it a dictionary"""

        arr = np.array(coll)
        dimN = arr.shape[0]
        new_dict = {}

        for i in xrange(dimN):
            new_dict.update({i: self.filter_zero_prob(coll[i])})

        return new_dict


class markov_seq(markov_chain):
    def __init__(self, order, train_data=None, transition_matrix=None):
        """Markov_seq is an object used to make a random sequence that mimics
           the base composition of the training data. If no training
           data is provided, it builds a sequence with uniform base
           probabilities. This object treats the addition of each
           nucleotide as a state transition; therefore, the Markov
           chain becomes represents a randomly generated sequence
           (MC.get_chain()).

        """
        
        markov_chain.__init__(self, order)

        if train_data is not None:
            self.coll = train_data.upper()
            self.transition_matrix = self.build_tr_matrix()
        # elif transition_matrix is not None:
            # stuff to directly use a transition matrix
            # generate train_data using the transition matrix
        else:
            raise NameError("noInputData")
        
        self.chain = self.gen_start_state()

    def gen_start_state(self):
        if self.order == 0:
            n = 1
        else:
            n = self.order

        i = np.random.randint(len(self.coll)-n-1)
        return self.coll[i:i+n]

    def build_tr_matrix(self):
        tr_dict = {}
        n = self.order
        
        if n == 0:
            cnts = Counter(self.coll)
            for state in cnts.keys():
                tr_dict.update({tuple(state): cnts})
        else:
            for i in xrange(0, len(self.coll)-n):
                state = tuple(self.coll[i:i+n])
                next_state = self.coll[i+n]
                if state in tr_dict:
                    cur = tr_dict.get(state).get(next_state, 0)  # nested dict call
                    tr_dict[state].update({next_state: cur+1})
                else:
                    tr_dict[state] = {next_state: 1}

        # # normalize the transition matrix to between [0, 1]
        for state, tr in tr_dict.iteritems():
            next_states = tr.keys()
            cnts = tr.values()
            sigma = sum(cnts)
            props = map(lambda c: float(c)/sigma, cnts)
            tr_dict.update({state: dict(zip(next_states, props))})

        return tr_dict

    def gen_next_state(self, nnext=None):
        if nnext is None:
            nnext = 1

        if self.order == 0:
            n = 1
        else:
            n = self.order
            
        for _ in xrange(nnext):
            state = tuple(self.chain[-n:])
            tr = self.make_transition(self.transition_matrix[state])
            self.chain = self.chain + tr

    def get_chain(self):
        return self.chain

    def __str__(self):
        return self.get_chain()


class markov_evolution(markov_chain):
    def __init__(self, inseq, transition_matrix):
        """The markov_evolution object handles the sequence evolution
           logic. It knows the starting sequence and can perform a
           single generation evolution using the transition
           matrix. The mutation method is defined by the mutate_seq()
           function, which returns a new sequence to represent the
           next generation.
        """

        markov_chain.__init__(self, order=0,
                              transition_matrix=transition_matrix)
        self.inseq = list(inseq)

    # Takes a sequence and mutates it on at an average rate of lambda = mu
    def mutate_seq(self, mu):
        """Takes a mean number of mutations per sequence (mu), mutation model,
        and a sequence, returns a new sequence that contains
        mutations. the number of mutations per sequence is determined
        by using a poisson distributed random number, typically on
        average there are (10**-9 * seq_length) nucleotides mutated.

        """
        
        mut_num = np.random.poisson(mu)
        n = len(self.inseq)

        # get the number of mutations to make and select the n positions
        # to mutate
        mut_pos = su.take(mut_num,
                          su.distinct(
                              su.repeatedly(np.random.randint, 0, n-1)))

        # applies the mutations to the sequence at each mutation position
        for i in mut_pos:
            orig_nt = self.inseq[i]
            new_nt = self.make_transition(self.transition_matrix[orig_nt])
            self.inseq[i] = new_nt

    def __str__(self):
        return ''.join(self.inseq)

    def as_str(self):
        return self.__str__()
