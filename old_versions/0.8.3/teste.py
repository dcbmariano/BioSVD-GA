
import numpy as np
from Bio import SeqIO
from random import randint
from BioSVD import *
from scipy import linalg

mat = np.zeros(shape=(2, 3))
mat[0][0] = 3
mat[0][1] = 2
mat[0][2] = 2

mat[1][0] = 2
mat[1][1] = 3
mat[1][2] = -1

#mat[2][0] = 0
#mat[2][1] = 1
#mat[2][2] = -1

print mat
u,s,v = linalg.svd(mat)

s = np.diag(s)

print "u",u, "\n"
print "s",s, "\n"
print "v",v.T, "\n"


'''

f = open("1000.fasta", "w")
seq_fasta = list(SeqIO.parse(open("arquivo2.fasta", "r"),"fasta"))

i=0
while i< 1000:
	t = randint(0,len(seq_fasta))
	SeqIO.write( seq_fasta[t], f , "fasta" )
	del seq_fasta[t]
	i= i+1
'''

