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
import os


# Testing function "read"
seqs = BioSVD.read('fasta',sys.argv[1:])

# Testing function "kmer"
mfkmer = BioSVD.kmer(seqs,3,'a')

# Testing function "svds"
[U, S, V] = BioSVD.svds(mfkmer)

# Testing function "factor"
s = BioSVD.factor(S,'plot')

# Testing function "extractFactor"
mfkmer_3 = BioSVD.extractFactor(S,V,3)

# Testeing function "reductor"
reduced_matrix = BioSVD.reductor(mfkmer_3 , 3)

# Testing function "plot2"
BioSVD.plot2(mfkmer_3,'plot')

# Testing function "plot3"
BioSVD.plot3(mfkmer_3,'plot')

# Distance all against all
dist = BioSVD.dist(mfkmer_3)

# Dalaunay
contacts = BioSVD.delaunay(mfkmer_3)

# Pairs of contacts
pairs = BioSVD.uniquePairs(contacts)

# Knearest, e.g.: for k = 3
nearest = BioSVD.knearest(dist, 1, 3)

print "The three nearest neighboors of %s are:\n1. %s\n2. %s\n3. %s\n" \
%(BioSVD.info(seqs,1,'name'), BioSVD.info(seqs,nearest[0],'name'), BioSVD.info(seqs,nearest[1],'name'), BioSVD.info(seqs,nearest[2],'name'))



