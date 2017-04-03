import numpy as np
from collections import Counter
import bisect


class markov_seq():
    def __init__(self, order, train_data = None, transition_matrix = None):
        """

        """
        self.order = order

        if train_data is not None:
            self.coll = train_data.upper()
            self.tr_dict = self.build_tr_matrix()
        elif transition_matrix is not None:
            # stuff to directly use a transition matrix
            # generate train_data using the transition matrix
        else:
            raise("provide either training data or a transition matrix")
        
        self.chain = self.start_state()

    def start_state(self):
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

    def transition(self, tr_arr):
        """Determine mutation based on transition probabilities"""

        r = np.random.rand()
        states = tr_arr.keys()
        cum_probs = np.cumsum(tr_arr.values())
        return states[bisect.bisect(cum_probs, r)]

    def gen_next_state(self, nnext=None):
        if nnext is None:
            nnext = 1

        if self.order == 0:
            n = 1
        else:
            n = self.order
            
        for _ in xrange(nnext):
            state = tuple(self.chain[-n:])
            tr = self.transition(self.tr_dict[state])
            self.chain = self.chain + tr

    def get_chain(self):
        return self.chain
