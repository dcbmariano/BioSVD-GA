#!/usr/bin/python
#     Program: biosvd.py
#    Function: Implementa o metodo de svd para classificacao de proteinas
# Description: inclui: delaunay.py e sort.py - Por enquanto estamos limitado a 23 cores no plot assim so podemos ter 23 familas diferentes. Isso sera revisado(Por Thiago)
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 1


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
nomePlotClusterizacao = ''
nomePlotPosto = ''

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

	if( sys.argv[i] == '-g'):
		nomePlotClusterizacao = sys.argv[i+1]

	if( sys.argv[i] == '-p'):
		nomePlotPosto = sys.argv[i+1]


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


A = np.zeros(shape=(10, 5))
A[0] = [0, 0, 0, 1, 0]
A[1] = [0, 0, 0, 0, 1]
A[2] = [0 ,0, 0, 0, 1]
A[3] = [1 ,0, 1, 0, 0]
A[4] = [1 ,0, 0, 0, 0]
A[5] = [0 ,1, 0, 0 ,0]
A[6] = [1 ,0, 1, 1, 0]
A[7] = [0 ,1, 1, 0, 0]
A[8] = [0 ,0, 1, 1 ,1]
A[9] = [0 ,1, 1, 0 ,0]

q = np.zeros(shape=(1, 10 ))
q[0]=[0, 0,0, 0, 0, 0, 0, 1, 1, 1]

#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
Cores = ['b', 'g', 'r', 'c','m','y', 'k', 'w', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
		'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

print "SVD"
aux,T = bio.SVD( A, 2 , nomePlotPosto )

print "Creating graphic 3D"
fig1 = plt.figure()
fig1.suptitle('Clusterizacao')

print "Model"
x = aux[ 0:1, :]
y = aux[ 1:2, :]
plt.scatter(x, y ,c=Cores[2], marker='o', s = 100, label = 'Model' )

print "Query"
T2 =  T[ :, 0:2]
qtil = np.dot(T2.T , q.T)
plt.scatter(qtil[0], qtil[1], c='#000000', marker='*',label='Query', s = 100 )
	
# Criando figura
plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=10, bbox_to_anchor=(0, 0))
if(len(nomePlotClusterizacao) > 0):
	fig1.savefig('results/'+nomePlotClusterizacao+'.png', dpi=300)
	plt.close(fig1)
else:
	plt.show()


# Fim do tempo de execucao 
fim = time.time()
print "Time: ", fim-ini
