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

os.system("clear")

# Controle do tempo de exucucao
ini = time.time()

# Parametros e help
i = 1

while i < (len(sys.argv)):

	# m = modelo
	if( sys.argv[i] == '-m'):
		FileModelo = sys.argv[i+1]

	# mt = modelo - arquivo tabular
	if( sys.argv[i] == '-mt'):
		FileModeloTabular = sys.argv[i+1]

	# q = query
	if( sys.argv[i] == '-q'):
		FileQuery = sys.argv[i+1]

	# qt = query - arquivo tabular
	if( sys.argv[i] == '-qt'):
		FileQueryTabular = sys.argv[i+1]

	# Help
	if sys.argv[i] == '-h' or sys.argv[i] == '--help':
		print "\n*** BIOSVD ***\nSyntax: \n\tpython biosvd.py \n\t\t-m  [model] \n\t\t-mt [model tabular file] \n\t\t-q  [query] \n\t\t-qt [query tabular file - optional]\n"
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


hashTabModelo, FamiliasModelo = bio.CriarHashTab(FileModeloTabular)
hashTabQuery, FamiliasQuery = bio.CriarHashTab(FileQueryTabular)

print "Parsing files"
sequenciasModelo = list(SeqIO.parse( open(FileModelo, "r")  , "fasta"))
sequenciasQuery = list(SeqIO.parse( open(FileQuery, "r")  , "fasta"))

print "Model number sequence: %d \nQuery number sequence: %d" %( len(sequenciasModelo), len(sequenciasQuery) )

print "Sorting Model Seq"
NumFamiliasModelo ,FamiliasModelo, DistribuicaFamiliasModelo , sequenciasModelo = bio.Sort( sequenciasModelo, hashTabModelo )

numeroSeqModel = len(sequenciasModelo)
numeroSeqQuery = len(sequenciasQuery)
NumroDeSequencias = numeroSeqModel + numeroSeqQuery


#Jutando todas as sequencias numa so lista
allSequences = []
for i in sequenciasModelo:
	allSequences.append(i)
for i in sequenciasQuery:
	allSequences.append(i)

#Preenchendo matriz de frequencia
matFrequencia = bio.Kmer( allSequences  ,3 )
print "Length of the matrix: ",matFrequencia.shape


#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
Cores = ['b', 'g', 'r', 'c','m','y', 'k', 'w', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
		'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

print "SVD"
aux = bio.SVD( matFrequencia, 3 )

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

	ax.scatter(x, y, z, c=Cores[IDCor], marker='o', label=FamiliasModelo[IDCor], s = 100 )
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
if(len(sys.argv) >= 10):
	fig1.savefig('results/'+str(sys.argv[10])+'.png', dpi=300)
	plt.close(fig1)
else:
	plt.show()

# Calculando delauney
print "Calculing delaunay"
bio.delaunay( temp, aux , sequenciasQuery, hashTabQuery)

print "| Running validation | %s"  %FileQueryTabular
bio.Validation( hashTabQuery, sequenciasQuery)

# Fim do tempo de execucao 
fim = time.time()
print "Time: ", fim-ini
