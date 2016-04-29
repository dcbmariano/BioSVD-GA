#!/usr/bin/python
# coding=UTF-8
#     Program: BioSVD.py
#    Function: Protein and nucleotide clusterization
# Description: Library with
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#        Date: April 18, 2016
#     Version: 0.8.1 alpha


# Requirements
from numpy import linalg as LA
from scipy.sparse.linalg import svds as SVDS
from Bio import SeqIO
import itertools
import numpy as np
from scipy.spatial import Delaunay,distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import os
from functools import partial
from random import *
import sys


# List of functions:
# - read: OK
# - kmer: OK
# - svds: OK
# - factor: OK
# - extractFactor: OK
# - reductor
# - plot2: OK
# - plot3: OK
# - model
# - query
# - knearest: OK
# - info: OK
# - voting
# - delaunay: OK
# - uniquePairs: OK
# - crossValidation


# Function read: allows reading of multiple files in FASTA or PDB format
# input: (1) format ('fasta' or 'pdb') - (2) list of files
def read(format,*files):
	
	seqs = []
	i = 0

	if format == 'fasta' or format == 'FASTA':
		for f_listas in files:
			for f in f_listas:
				seq_fasta = list(SeqIO.parse(open(f, "r"),"fasta"))

				for s in seq_fasta:
					seqs.append({'id':i, 'name': s.name, 'seq': str(s.seq), 'group': f})
					i += 1

	return seqs


# Function kmer: receives a sequence object and return a k-mer matrix
def kmer(seqs,kmer_len,type):

	# Creating K-MERS

	# For amino acids
	if type == 'a' or type == 'A':
		aminoacidos = [ 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]
		HashKMer = {}
		t = 0
		kmers = [''.join(p) for p in itertools.product(aminoacidos, repeat=kmer_len)]

		NumerosDeKMER = len(kmers)
		numseqs = len(seqs)

		matFreq = np.zeros(shape=(NumerosDeKMER, numseqs))
		
		# Preenchendo a tabela HASH com os Kmers
		for i in kmers:
			HashKMer[ i ] = t
			t = t+1

		# Filling K-MER Matrix
		for j in range(numseqs):
			for k in range( len(seqs[j]['seq']) -(kmer_len-1)):
				#Kmer contem algum aminoacidos fora dos 20 que estamos usando. Vamos ignorar tal kmer
				if seqs[j]['seq'][k:k+kmer_len] in HashKMer: # O preenchimento e por coluna
					i = HashKMer.get(seqs[j]['seq'][k:k+kmer_len])
					matFreq[i][j] = 1

		return matFreq

	# For nucleotides
	elif type == 'n' or type == 'N':
		nucleotides = ['A','T','C','G']
		HashKMer = {}
		t = 0
		kmers = [''.join(p) for p in itertools.product(nucleotides, repeat=kmer_len)]

		NumerosDeKMER = len(kmers)
		numseqs = len(seqs)

		matFreq = np.zeros(shape=(NumerosDeKMER, numseqs))
		
		# Preenchendo a tabela HASH com os Kmers
		for i in kmers:
			HashKMer[ i ] = t
			t = t+1

		# Filling K-MER Matrix
		for j in range(numseqs):
			for k in range( len(seqs[j]['seq']) -(kmer_len-1)):
				if seqs[j]['seq'][k:k+kmer_len] in HashKMer: 
					i = HashKMer.get(seqs[j]['seq'][k:k+kmer_len])
					matFreq[i][j] = 1

		return matFreq


def svds(matrix):
	k = matrix.shape[1]
	U ,s, Vt = SVDS(matrix, k-1)
	S = np.diag(s)	# o svds retorna um vetor com os valores singulares da matriz. Aqui criamos uma matriz diagonal
					# com os valores singulares. O Script exige que tenahmos uma matriz  e não um vetor.
	V = Vt.transpose()	#O svds retorna uma matriz V que e a transposta da matriz retornada pelo matlab. Entao aqui 						#corrigimos isso
	return U,S,V


def factor(S,action):

	# Actions:
	# - plot: shows the factor figure  
	# - save: save the factor figure in the disk

	s = np.diag(S)
	s = s*(1/np.sum(s))
	s = list(reversed(s))

	# Generating figure
	fig = plt.figure()
	fig.suptitle('Factor plot')
	plt.plot(s)

	# Save figure
	if action == 'save':
		fig.savefig('factor.png', dpi=300)
		plt.close(fig)

	# Ploting figure
	elif action == 'plot':
		plt.show()


def extractFactor(S,V,factor):
	Sfactor = S[:factor,:factor]
	Vfactor = V[:,:factor] 
	AUXfactor = np.dot(Sfactor, Vfactor.transpose()) 
	return AUXfactor


def reductor(matrix,dim):
	print "Coming soon."


def plot2(matrix,action):
	# Actions:
	# - plot: shows the factor figure  
	# - save: save the factor figure in the disk
	colors = np.random.rand(matrix.shape[1])
	s = matrix[:2,:]

	# Generating figure
	fig = plt.figure()
	fig.suptitle('Factor plot')
	x = s[:1,:]
	y = s[1:2,:]
	plt.scatter(x, y, s = 100, c=colors, alpha=0.5)

	# Save figure
	if action == 'save':
		fig.savefig('factor2d.png', dpi=300)
		plt.close(fig)

	# Ploting figure
	elif action == 'plot':
		plt.show()


def plot3(matrix,action):
	# Actions:
	# - plot: shows the factor figure  
	# - save: save the factor figure in the disk
	colors = np.random.rand(matrix.shape[1])
	s = matrix[:3,:]

	# Generating figure
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	fig.suptitle('Factor plot')
	x = s[:1,:]
	y = s[1:2,:]
	z = s[2:3,:]
	ax.scatter(x, y, z, c=colors, marker='o', s = 100)


	# Save figure
	if action == 'save':
		fig.savefig('factor3d.png', dpi=300)
		plt.close(fig)

	# Ploting figure
	elif action == 'plot':
		plt.show()


def model():
	print "Coming soon."


def query():
	print "Coming soon."


def knearest(dist, element, k):
	matrix = np.argsort(dist[element])

	return matrix [1:k+1]


def voting():
	print "a"


def info(seqs, position, info):
	for s in seqs:
		if s['id'] == position:
			return s[info]


def delaunay(matrix):

	element = []
	coord = []

	for i in range(len(matrix[0])):
		for j in range(len(matrix)):
			element.append(float(matrix[j][i]))

		coord.append(element)
		element = list()

	delaunay = Delaunay(coord)
	return delaunay.vertices


def uniquePairs(delaunay):
	pairs = []
	dimensions = len(delaunay[0])

	for pair in delaunay: 
		for k in range(dimensions + 1):
			for j in range(k+1,dimensions):
				pairs.append([pair[k],pair[j]])

	# Remove duplications ********************************************
	unique_pairs = []
	for x in pairs:
		if x not in unique_pairs:
			unique_pairs.append(x)

	return unique_pairs


def crossValidation():
	print "Coming soon."


def dist(matrix):
	# input matrix => m (rows) x n (cols)
	# output matrix => n x n (distance all against all)

	from math import sqrt

	elements = len(matrix[0])
	dist = []
	last = []
	d = 0

	for i in range(elements):
		for j in range(elements):
			# euclidian distance
			d = sqrt((matrix[0][i] - matrix[0][j])**2 + (matrix[1][i] - matrix[1][j])**2 + (matrix[2][i] - matrix[2][j])**2)
			last.append(d)

		dist.append(last)
		last = list()

	return dist
