# Seq-sim

## Introduction
The aim of Seq-sim is to help further our understanding of the
selective pressure to maintain a protein coding sequences. The central
dogma of biology is that DNA -> RNA -> Protein. Therefore, given a
starting protein coding sequence and a mutation rate, seq-sim
simulates mutations accumulating over many generations. If at any
point, the protein coding portion is disrupted, the sequence is
considered "dead" or non-functional. If the given sequence is known to
be functional, the ratio of functional to non-functional sequences
measures the selective pressure on the sequence. A ratio > 1 indicates
low selective pressure (genetic drift); a ratio < 1 indicates pressure
to maintain the coding sequence.

Seq-sim simulates sequence evolution using a Markov process. The user
needs to specify the mutation rate and the number of generation to
simulate. The Markov process is At each generation, we use a Poisson
random number to determine how many mutations to apply. For each
mutation, we assume that the base can mutate to any other base with
equal probability. The output of the program is the number of
funtional and non-functional sequences.

## How to run the program
python seq_generator.py

or

python seq_generator.py "ATGATCATACATGACAGGCTGCTTGGCGAATTCTACGTCAGTACACACCAAGGCTCTGCGCCCGCTGTCGAAAGCGCCTATCGCTAATGTCTGCTGTGGCGCATT"

