#!/usr/bin/python
#     Program: example.py
#    Function: Test library BioSVD
# Description: It shows tests using BioSVD for sequence classification
#      Author: TS Correia, DCB Mariano, JR Barroso, RC de Melo-Minardi
#     Version: 0.8.1


# Help - Use: python example.py file_1.fasta file_2.fasta file_n.fasta


# Requirements
from BioSVD import *
import sys


# Testing function "read"
seqs = BioSVD.read('fasta',sys.argv[1:])

# Testing function "kmer"
mfkmer = BioSVD.kmer(seqs,3,'a')

# Testing function "svd"
[U, S, V] = BioSVD.svd(mfkmer)

# Testing function "factor"
s = BioSVD.factor(S,'plot')

# Testing function "extractFactor"
mfkmer_3 = BioSVD.extractFactor(S,V,3)