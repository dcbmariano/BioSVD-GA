#!/usr/bin/python
#     Program: tolerancia.py
#    Function: Test library BioSVD
# Description: It shows tests using BioSVD for sequence classification
#      Author: TS Correia, DCB Mariano, JR Barroso, RC de Melo-Minardi
#     Version: 0.8.3 alpha

# Help - Use: python tolerancia.py file_1.fasta file_2.fasta file_n.fasta


# Requirements
from BioSVD import *
import numpy as np
import sys
import os
os.system("clear")

#Numero de tolerantes conhecidas
n = 21

# Reading the files
print "Reading"
seqs = BioSVD.read('fasta',sys.argv[1:])
print "Number of sequences", len(seqs)
# Creating kmer
print "kmer"
mfkmer = BioSVD.kmer(seqs,3,'a')

# Applying SVDS
print "SVD"
[U, S, V] = BioSVD.svds(mfkmer)

#Printing Factor
#factor = BioSVD.factor(S,"save")
factor = 3
print "Factor: ", factor

# "extractFactor"
print "extractFactor"
mfkmer_factor = BioSVD.extractFactor(S,V,factor)

#Reduzindo a dimensao para 3
#mat_reduced = BioSVD.reductor(mfkmer_factor.T ,3)

print "Calculating centroid"
coord,r = BioSVD.centroide(mfkmer_factor[: , :n],factor )

# Plotting
BioSVD.plot2(mfkmer_factor,'save', n ,coord,r  )
BioSVD.plot3(mfkmer_factor,'save', n )


print "Radius:", r
seqPosition = BioSVD.centroideSeqs(coord,r, mfkmer_factor, n)
print "Number of predicted sequences",len(seqPosition)

print "Saving predicted sequences"
BioSVD.save(seqs, seqPosition)























