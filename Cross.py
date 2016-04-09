#!/usr/bin/python
#     Program: biosvd.py
#    Function: Implementa o metodo de svd para classificacao de proteinas
# Description: inclui: delaunay.py e sort.py - Por enquanto estamos limitado a 23 cores no plot assim so podemos ter 23 familas diferentes. Isso sera revisado(Por Thiago)
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 7.01


# Hierarquia do programa ******************************************
#
# ./        	Diretorio raiz
# biosvd.py 	Script principal
# README.txt	Contem instruoes de execucao
# example   	Diretorio onde ficaram armazeados dados para testes
# results   	Diretorio onde deverao ser salvos os outputs	 


# IMPORTS
from Biosvd import biosvd as bio
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import os
from random import *

os.system("clear")
try:
	os.mkdir("results")
except:
	print "Erro ao criar diretorii results"

# Controle do tempo de exucucao
ini = time.time()

# Parametros e help
i = 2
nomePlotClusterizacao = ''
nomePlotPosto = ''
kfold = 0

opcao_entrada =  sys.argv[1]

while opcao_entrada == '-B' and i< 8:
	if( sys.argv[i] == '-g'):
		nomePlotClusterizacao = sys.argv[i+1]

	if( sys.argv[i] == '-p'):
		nomePlotPosto = sys.argv[i+1]

	if( sys.argv[i] == '-k'):
		kfold = int(sys.argv[i+1])
	i = i+2


while opcao_entrada == '-A' and i < (len(sys.argv)):

	# m = modelo
	if( sys.argv[i] == '-m'):
		FileSequencias = sys.argv[i+1]

	# mt = modelo - arquivo tabular
	if( sys.argv[i] == '-mt'):
		FileTabular = sys.argv[i+1]

	# k-fold
	if( sys.argv[i] == '-k'):
		kfold = int(sys.argv[i+1])

	if( sys.argv[i] == '-g'):
		nomePlotClusterizacao = sys.argv[i+1]

	if( sys.argv[i] == '-p'):
		nomePlotPosto = sys.argv[i+1]


	# Help
	if sys.argv[i] == '-h' or sys.argv[i] == '--help':
		print "\n*** BIOSVD ***\nSyntax: \n\tpython biosvd.py \n\t\t -A \n\t\t-m  [model] \n\t\t-mt [model tabular file] \n\t\t-q  [query] \n\t\t-k [K-fold]\n"
		sys.exit()
	
	i = i+2

# Valida se nenhum parametro for enviado
if i == 1:
	print "\n*** BIOSVD ***\nSyntax: \n\tpython biosvd.py \n\t\t-m  [model] \n\t\t-mt [model tabular file] \n\t\t-q  [query] \n\t\t-qt [query tabular file - optional]\n"
	sys.exit()


# Aqui comeca a magia
print "\n*******************************"
print "        *** BIOSVD ***"
print "*******************************\n"

hashTab = {}
allSequences = []
print "Parsing files"
#Criando tabular
if opcao_entrada == '-A':
	hashTab, FamiliasModelo = bio.CriarHashTab(FileTabular)
	allSequences = list(SeqIO.parse( FileSequencias , "fasta"))
else:
	for j in range(i,len(sys.argv)):
		seqtemp = list(SeqIO.parse( sys.argv[j], "fasta"))
		allSequences = allSequences + seqtemp
		for p in seqtemp:
			hashTab[ p.id ] = sys.argv[j].split('.')[0]


print bio.CrossValidation( allSequences,kfold, hashTab ,opcao_entrada)
'''
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

	NumFamiliasModelo ,FamiliasModelo, DistribuicaFamiliasModelo , sequenciasModelo =  bio.Sort( sequenciasModelo, hashTab ,opcao_entrada )

	print "Model lenth: %d\nQuery lenth %d:\n" % (numeroSeqModel, numeroSeqQuery )


	#Preenchendo matriz de frequencia
	matrizFrequencia = bio.Kmer( sequenciasModelo +sequenciasQuery  , 3 )
	#Q = bio.Kmer( sequenciasQuery  , 3 )
	print "Length of the matrix Model: ",matrizFrequencia.shape
	#print "Length of the matrix Query: ",Q.shape

	#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
	#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
	Cores = ['b', 'g', 'r', 'c','m','y', 'k', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
			'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

	print "SVD"
	aux ,T= bio.SVD( matrizFrequencia, 3 , nomePlotPosto )

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
	# POG DIEGO - modificar isso no futuro *************************************************************************
	print "Calculing delaunay"

	bio.delaunay( temp , aux, sequenciasQuery,hashTab, opcao_entrada )
	print "| Running validation |"
	A,B = bio.Validation(hashTab, sequenciasQuery, opcao_entrada )
	acertos = acertos +  A
	erros = erros +B

print "Acuracia: %0.3f" % float(acertos/float((erros+ acertos)))
'''
# Fim do tempo de execucao 
fim = time.time()
print "Time: ", fim-ini
