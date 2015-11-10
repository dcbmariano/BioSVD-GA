#!/usr/bin/env python
# -*- coding: latin-1 -*-


from numpy import linalg as LA
from Bio import SeqIO
import itertools
import numpy as np
from scipy.spatial import Delaunay,distance

def SVD (  matriz ):

	U, s, V = LA.svd( matriz )
	return U, s, V.transpose()
 
def Kmer( sequencias, k ):
	print "Creating K-MERS"
	aminoacidos = [ 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]
	HashKMer = {}
	t = 0
	k=3
	kmers = [''.join(p) for p in itertools.product(aminoacidos, repeat=k)]

	NumerosDeKMER = len(kmers)
	NumroDeSequencias = len(sequencias)

	matFrequencia = np.zeros(shape=(NumerosDeKMER, NumroDeSequencias))
	
	# Preenchendo a tabela HASH com os Kmers
	for i in kmers:
		HashKMer[ i ] = t
		t = t+1

	print 'Filling K-MER Matrix'
	for j in range(NumroDeSequencias):
		for k in range( len(sequencias[j]) -2):
			#Kmer contem algum aminoacidos fora dos 20 que estamos usando. Vamos ignorar tal kmer
			if str(sequencias[j].seq[k:k+3]) in HashKMer: # O preenchimento e por coluna
				i = HashKMer.get(str(sequencias[j].seq[k:k+3]))
				matFrequencia[i][j] = 1

	return matFrequencia


#Retorna onde ocorre a maior variacao dos elementos da diagonal da matriz S
def Posto( matrizS ):
	p = 1
	S =  np.diag(matrizS)

	#Normlizando a matriz S
	Temp = S/np.linalg.norm(S)
	SN = np.diag( Temp ) # matriz normalizada

	
	
	Delta = SN[1] - SN[0]
	
	if Delta < 0:
		Delta = -1* Delta
	#print Delta

	for i in range( len(SN) -1 ):
		temp = SN[i+1] - SN[i]
		if temp < 0:
			temp = -1* temp

		#print temp

		if temp > Delta:
			p = i+1
			Delta = temp

	return p
