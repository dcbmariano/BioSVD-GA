#!/usr/bin/python
#     Program: example.py
#    Function: Test library BioSVD
# Description: It shows tests using BioSVD for sequence classification
#      Author: TS Correia, DCB Mariano, JR Barroso, RC de Melo-Minardi
#     Version: 0.8.3 alpha



# Requirements
from BioSVD import *
import numpy as np
import sys
import os
os.system("clear")


# Reading the files
print "Reading"
mat = BioSVD.read('csv',sys.argv[1:])
mat = mat.transpose()
print "Number of sequences", mat.shape

# Applying SVDS
print "SVD"
[U, S, V] = BioSVD.SVD(mat)

#Printing Factor
BioSVD.factor(S,"plot", "svd")
factor = int(raw_input('Entre com o factor: '))

# "extractFactor"
print "extractFactor"
mat_factor = BioSVD.extractFactor(np.diag(S),V,factor)

BioSVD.saveCSV(mat_factor, "out.csv")
























