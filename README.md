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
simulate. At each generation, we use a Poisson random number to
determine how many mutations to apply. For each mutation, apply the
Jukes-Cantor evolutionary model, which assumes that the base can
mutate to any other base with equal probability. Other models such as
HKY85 and GTR will be implemented in future versions.

The final output of the program is the number of funtional and
non-functional sequences. Functional in the current version is defined
as a sequence without premature stop codons.

## How to run the program
```bash
python seq_generator.py -h
```

or

```bash
python seq_generator.py -i "ATGATCATACATGACAGGCTGCTTGGCGAATTCTACGTCAGTACACACCAAGGCTCTGCGCCCGCTGTCGAAAGCGCCTATCGCTAATGTCTGCTGTGGCGCATT"
```

`Counter({True: 93, False: 7})`
This result implies that a sequence of the given sequence simulated through 1000 generations, and repeated 100 times, will be functional (no premature stops) 93% of the time. This can be effectively thought of as answering the question, "given a starting sequence, repeat evolution 100x, how often will the gene remain functional?"

## TODO
- Add other evolutioary models (HKY85, GTR)
- Simulate indels
- Write more test cases
- Clean up main function in seq_generator and stop printing test cases