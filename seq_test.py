import unittest
import seq_generator
import seq_markov
import seq_object

nmem = 1
inseq = "CATACATGACAGGCTGCTTGGCGAATTCTACGTCAGTACACACCAAGGCTCTGCGCCCGCTGTCGAAAGCGCCTATCGCTAATGTCTGCTGTGGCGCATT"


class seq_tests(unittest.TestCase):

    def test_chain(self):
        MC = seq_markov.markov_seq(nmem, inseq)
        for _ in xrange(10):
            MC.gen_next_state()
        self.assertEqual(len(MC.get_chain()), 11)

    def test_tr_dict(self):
        MC = seq_markov.markov_seq(1, str('A')*100)
        self.assertEqual(MC.tr_dict, {('A',): {'A': 1.0}})

    def test_tr_dict2(self):
        MC = seq_markov.markov_seq(0, 'ACTGACTG')
        self.assertEqual(MC.tr_dict[tuple('A')],
                         {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25})

    def test_gen_next(self):
        MC = seq_markov.markov_seq(0, str('A')*100)
        MC.gen_next_state()
        updated_chain = MC.get_chain()
        self.assertEqual('AA', updated_chain)

    def test_gen_next2(self):
        MC = seq_markov.markov_seq(1, 'ACGTACGT')
        MC.chain = 'A'
        MC.gen_next_state(11)
        self.assertEqual('ACGTACGTACGT', MC.get_chain())


if __name__ == "__main__":
    unittest.main()
