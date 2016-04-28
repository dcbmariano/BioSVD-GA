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
# - plot2
# - plot3
# - model
# - query
# - knearest
# - kmeans
# - delaunay
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
	print s
	s = s*(1/np.sum(s))

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
	Vfactor = V[:,:factor] #Modifiquei aqui pois na apostila está assim 
	AUXfactor = np.dot(Sfactor, Vfactor.transpose()) #Modifiquei aqui pois na apostila está assim 
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


def knearest():
	print "Coming soon."


def kmeans():
	print "Coming soon."


def delaunay(matrix, action):

	# Actions:
	# - invert: invert columns and lines in the matrix

	if action == 'invert':
		print "Coming soon."

	return Delaunay(matrix)


def crossValidation():
	print "Coming soon."


def dist():
	print "Coming soon."


def cosine(A,B):
	dot = A*B 
	normA = LA.norm(A)
	normB = LA.norm(B)

	cos = float( dot / ( normA * normB) )
	return cos


# ---------------------------------------------------------------------------------------------------------------------


#Calculates euclidean distances
def dist_euclidiana(p1,p2):
	#print "dist_euclidiana"
	#print p1.shape
	#print p2.shape 
	return LA.norm(p1-p2)

#Calculamos as 5 mais proximas proteinas em relacao a j-esima proteina
def MaisProxima( matPositionModel , numeroSeqModel,  Query ):
	ID = 0
	aux = 0
	dist = 0
	
	aux = dist_euclidiana( Query , matPositionModel[: ,0:1] )
	for i in range( 1, numeroSeqModel):
		#print i
		dist = dist_euclidiana( Query , matPositionModel[:, i:i+1] )
		if dist < aux:
			aux = dist
			ID = i
	return ID

def DistanciaEuclidina(matPositionModel ,matPositionQuery, NumroDeSequenciasQuery, numeroSeqModel, sequenciasModelo ,sequenciasQuery):

	family_predicted = {}
	'''
	pool = mp.Pool()
	func = partial(MaisProxima, matPositionModel, numeroSeqModel )

	IDS = pool.map(func, matPositionQuery.T )
	print "IDS"
	print len(IDS)
	'''

	for i in range( NumroDeSequenciasQuery ):
		ID = MaisProxima( matPositionModel, numeroSeqModel,  matPositionQuery[:, i:i+1] )
		IDQuery = sequenciasQuery[i].id.split("|")[2]
		IDpredicted = sequenciasModelo[ID].id.split("|")[2]
		family_predicted[ IDQuery  ] = IDpredicted #ID da query / ID da familia predicted

	return family_predicted


def CrossValidation( allSequences, kfold ,hashTab ,opcao_entrada):

	acertos = 0
	erros =0
	numeroSequencias = len(allSequences)
	numero_Seq_por_grupos= numeroSequencias/kfold
	numeroSeqModel = numero_Seq_por_grupos
	listSequencias = [] #Contem o Kfold grupos de sequencias. cada grupo tem o mesmo numero de sequencias

	print "Creating K-fold"
	for j in range( kfold ):
		seq = []
		while len(seq) < numero_Seq_por_grupos :
			try:
				i =randint(0, len(allSequences))
				seq.append(allSequences[i])
				del allSequences[i]
			except:
				PO = 0
		listSequencias.append(seq[:])
		del seq[:]

	l=0
	Allgrups = listSequencias[:]
	print '\n'
	for sequenciasModelo in Allgrups:
		print "Grupo: %d" % (l+1)
		nomePlotPosto ="posto%d" % l
		nomePlotClusterizacao="Clus%d" % l
		l = l +1

		listSequencias.remove(sequenciasModelo) #Removendo o grupo de sequencias que sera o modelo
		sequenciasQuery = []
		for aux in listSequencias:
			sequenciasQuery = sequenciasQuery + aux

		numeroSeqQuery = len(sequenciasQuery)
		NumroDeSequencias = numeroSeqModel + numeroSeqQuery
		listSequencias.append(sequenciasModelo)
		print "Sorting Seq"
		# Numero de familas do modelo | Lista com todas as familias sem repeticao | Lista com a distruicao das Familias na lista | Lista com as sequencias do modelo

		NumFamiliasModelo ,FamiliasModelo, DistribuicaFamiliasModelo , sequenciasModelo =  Sort( sequenciasModelo, hashTab ,opcao_entrada )

		print "Model lenth: %d\nQuery lenth %d:\n" % (numeroSeqModel, numeroSeqQuery )

		#Preenchendo matriz de frequencia
		matrizFrequencia = Kmer( sequenciasModelo +sequenciasQuery  , 3 )
		print "Length of the matrix Model: ",matrizFrequencia.shape

		#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
		#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
		Cores = ['b', 'g', 'r', 'c','m','y', 'k', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
				'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

		print "SVD"
		aux ,T= SVD( matrizFrequencia, 3 , nomePlotPosto )

		print "Creating graphic 3D"
		fig1 = plt.figure()
		ax = fig1.add_subplot(111, projection='3d')

		F = {}
		tF = {}
		IDCor = 0
		temp = []

		print "Model"
		for i in range(len(DistribuicaFamiliasModelo)):
			if i == 0:#Primeira Familia Modelo
				x = aux[0:1,0:int(DistribuicaFamiliasModelo[0]) ]
				y = aux[1:2,0:int(DistribuicaFamiliasModelo[0]) ]
				z = aux[2:3,0:int(DistribuicaFamiliasModelo[0]) ]
				F[FamiliasModelo[i]] = "0:"+str(int(DistribuicaFamiliasModelo[0]) )
		
			elif(i == len(DistribuicaFamiliasModelo)-1 ):#Ultima familia Modelo
	
				x = aux[0:1,int(DistribuicaFamiliasModelo[ -1]) : numeroSeqModel ]
				y = aux[1:2,int(DistribuicaFamiliasModelo[ -1]) : numeroSeqModel ]
				z = aux[2:3,int(DistribuicaFamiliasModelo[ -1]) : numeroSeqModel ]
				F[FamiliasModelo[i]] = str(int(DistribuicaFamiliasModelo[i-1]) ) + ":" + str(numeroSeqModel)
		
			else:
				x = aux[0:1,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
				y = aux[1:2,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
				z = aux[2:3,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
				F[FamiliasModelo[i]] = str(int(DistribuicaFamiliasModelo[i-1]))+":"+str(int(DistribuicaFamiliasModelo[i]))

			#ax.scatter(x, y, z, c=Cores[IDCor], marker='o', label=FamiliasModelo[IDCor], s = 100 )
			IDCor = IDCor +1
			temp.append(FamiliasModelo[i])
			temp.append(str(F[FamiliasModelo[i]]))


		print "Query"
		tx = aux[0:1, numeroSeqModel +1: ]
		ty = aux[1:2, numeroSeqModel +1: ]
		tz = aux[2:3, numeroSeqModel +1: ]
		tF['Query family'] = str(numeroSeqModel+1)+":"+str(NumroDeSequencias )
		ax.scatter(tx, ty, tz, c='#000000', marker='*',label='Query', s = 100 )

		# Criando figura
		plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=10, bbox_to_anchor=(0, 0))
		if(len(nomePlotClusterizacao) > 0):
			fig1.savefig('results/'+nomePlotClusterizacao+'.png', dpi=300)
			plt.close(fig1)
		else:
			plt.show()

		# Calculando delauney
		print "Calculing delaunay"
		delaunay( temp , aux, sequenciasQuery,hashTab, opcao_entrada )

		print "| Running validation |"
		A,B = Validation(hashTab, sequenciasQuery, opcao_entrada )
		acertos = acertos +  A
		erros = erros +B

	return "%0.3f" % float(acertos/float((erros+ acertos)))










