# paHMM-dist
paHMM-Tree - pairwise statistical phylogenetic distance estimation using pair hidden Markov models. 

Forked from Marcin Bogusz' and Simon Whelan's paHMM-Tree project.

paHMM-Tree, pronounced as 'palm tree' is a phylogenetic tree inference tool that uses the principle of pairwise statistical alignment to construct a neighbour joining trees from raw sequences.

What paHMM-Tree does

paHMM-Tree outputs phylogenetic trees in NEWICK format from raw sequence data in FASTA format. It does not require a prior sequence alignment. If you feed it with aligned sequenced, the gaps will be ignored. Being based on a pair hidden Markov model (pair-HMM) paHMM-Tree is able to integrate over all possible alignments between all the pairs of sequences in the data set. By doing so it accounts for alignment uncertainty.

What it does not do

paHMM-Tree does not align the sequences. In principle it could, but the whole idea is to avoid conditioning on a single multiple sequence alignment (MSA). Also constructing MSAs from pairwise alignments might not lead to reliable results.

Code and compilation

The code is available on GitHub repository and accessible using the link above. The sources come with an Eclipse CDT project and a Makefile. The makefile should be good for most of the modern x86-based architectures. If your system does not support SSE2, change the compiler flags.

Documentation and Binaries

Check the project website at pahmm-tree.tk for the latest documentation and binaries 
