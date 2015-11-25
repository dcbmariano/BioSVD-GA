#!/usr/bin/python
#
#     Program: Experimento.py
#    Function: Implementa o metodo de svd para classificacao de proteinas
# Description: inclui: delaunay.py e sort.py
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 1
#OBS: Por enquanto estamos limitado a 23 cores no plot assim so podemos ter 23 familas diferentes. Isso sera revisado(Por Thiago)

# Generalizando
#python Experimento.py -m Experimento/1000_v1.fasta -mt Experimento/tabular_v1.tab -g plotClus


# IMPORTS
from Biosvd import biosvd as bs

from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg as LA
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import sys
import os
from os import mkdir
from random import *

os.system("clear")
try:
	mkdir("results")
except:
	print ""




# Controle do tempo de exucucao
ini = time.time()

# Declaracoes iniciais
FileSequencias = ''
FileSequenciasTab = ''
sequenciasQuery = []
#FileNomeQuery = ''
#FileNomeFamiliaQuery = ''


# Parametros e help
i = 1

while i < (len(sys.argv)):

	# m = modelo
	if( sys.argv[i] == '-m'):
		FileSequencias = sys.argv[i+1]

	# mt = modelo - arquivo tabular
	if( sys.argv[i] == '-mt'):
		FileSequenciasTab = sys.argv[i+1]

	# Help
	if sys.argv[i] == '-h' or sys.argv[i] == '--help':
		print "\n*** BIOSVD ***\nSyntax: \n\tpython Experimento.py \n\t\t-m  [model] \n\t\t-mt [model tabular file]\n"
		sys.exit()
	
	i = i+2

# Valida se nenhum parametro for enviado
if i == 1:
	print "\n*** BIOSVD ***\nSyntax: \n\tpython Experimento.py \n\t\t-m  [model] \n\t\t-mt [model tabular file]\n"
	sys.exit()


# Aqui comeca a magia
print "\n*******************************"
print "        *** BIOSVD ***"
print "*******************************\n"

print "Parsing Seq"
sequenciasModelo = list(SeqIO.parse( FileSequencias, "fasta"))

#temporario
#sequenciasModelo = list(SeqIO.parse( open("./Experimento/modelo.fasta", "r")  , "fasta"))
#sequenciasQuery = list(SeqIO.parse( open("./Experimento/query.fasta", "r")  , "fasta"))

HashTabular, Familias = bs.CriarHashTab(FileSequenciasTab)

# Clusterizacao 
print "Starting clusterization"	
# Declaracoes iniciais de clusterizacao 


#Verificar
# Selecionando a query
while len(sequenciasQuery) <= 264:
	try:
		i =randint(0, len(sequenciasModelo))
		sequenciasQuery.append( sequenciasModelo[i] )
		del sequenciasModelo[i]
	except:
		print ""


tamanhosequenciasModelo = len(sequenciasModelo)
tamanhosequenciasQuery = len(sequenciasQuery)
NumroDeSequencias = tamanhosequenciasModelo + tamanhosequenciasQuery

print "Sorting Seq"
# Numero de familas do modelo | Lista com todas as familias sem repeticao | Lista com a distruicao das Familias na lista | Lista com as sequencias do modelo
NumFamiliasModelo ,FamiliasModelo, DistribuicaFamiliasModelo , sequenciasModelo =  bs.Sort(sequenciasModelo, HashTabular )


print "Model lenth: %d\nQuery lenth %d:\n" % (tamanhosequenciasModelo, tamanhosequenciasQuery )


print 'Joining lists'
#Jutando todas as sequencias numa so lista
allSequences =  sequenciasModelo + sequenciasQuery


SeqIO.write(sequenciasQuery, open("./results/query.fasta", "w"), "fasta")
SeqIO.write(sequenciasModelo, open("./results/modelo.fasta", "w"), "fasta")


#Preenchendo matriz de frequencia
matrizFrequencia = bs.Kmer( allSequences  , 3 )
print "Length of the matrix: ",matrizFrequencia.shape

#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
Cores = ['b', 'g', 'r', 'c','m','y', 'k', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
		'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

print "SVD"
U, s, V = bs.SVD( matrizFrequencia )

K = 3
SK =  np.diag(s)

#Normlizando a matriz S
Temp = SK/np.linalg.norm(SK)
SN = np.diag( Temp )


fig1 = plt.figure(1)
fig2 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')

#Aqui limito os eixos X e Y do plot
plt.axis([0, 200, 0, 1])
plt.plot(SN.T)
fig2.savefig('./results/posto.png', dpi=300)
plt.close(fig2)

SK = SK[0:K,0:K]
VK = V[:,:K]
aux = np.dot(SK , VK.transpose() )

np.savetxt("MatrizPosicoes.txt", aux,fmt = '%.5f' ,delimiter = ';')


print "Creating graphic 3D"

# Salvando familias do modelo
#model = open("tmp/model.txt","w")
F = {}
tF = {}
IDCor = 0
temp = []
FamiliaCor = {}
k=0
x = []
y = []
z = []

print "Model"

for i in range(len(DistribuicaFamiliasModelo)):
	if i == 0:#Primeira Familia Modelo
		x = aux[0:1,0:int(DistribuicaFamiliasModelo[0]) ]
		y = aux[1:2,0:int(DistribuicaFamiliasModelo[0]) ]
		z = aux[2:3,0:int(DistribuicaFamiliasModelo[0]) ]
		F[FamiliasModelo[i]] = "0:"+str(int(DistribuicaFamiliasModelo[0]) )
		
	elif(i == len(DistribuicaFamiliasModelo)-1 ):#Ultima familia Modelo
	
		x = aux[0:1,int(DistribuicaFamiliasModelo[ -1]) : tamanhosequenciasModelo ]
		y = aux[1:2,int(DistribuicaFamiliasModelo[ -1]) : tamanhosequenciasModelo ]
		z = aux[2:3,int(DistribuicaFamiliasModelo[ -1]) : tamanhosequenciasModelo ]
		F[FamiliasModelo[i]] = str(int(DistribuicaFamiliasModelo[i-1]) ) + ":" + str(tamanhosequenciasModelo)
		
	else:
		x = aux[0:1,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
		y = aux[1:2,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
		z = aux[2:3,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
		F[FamiliasModelo[i]] = str(int(DistribuicaFamiliasModelo[i-1]))+":"+str(int(DistribuicaFamiliasModelo[i]))

	ax.scatter(x, y, z, c=Cores[IDCor], marker='o', label=FamiliasModelo[IDCor].strip() )
	IDCor = IDCor +1
	#print F[FamiliasModelo[i]]
	#familia = "%s\n%s\n" %(FamiliasModelo[i],F[FamiliasModelo[i]])
	#model.write(familia)
	temp.append(FamiliasModelo[i])
	temp.append(str(F[FamiliasModelo[i]]))

print "Query"
tx = aux[0:1, tamanhosequenciasModelo -1: ]
ty = aux[1:2, tamanhosequenciasModelo -1: ]
tz = aux[2:3, tamanhosequenciasModelo -1: ]
tF['Query family'] = str(tamanhosequenciasModelo-1)+":"+str(NumroDeSequencias )
ax.scatter(tx, ty, tz, c='#000000', marker='*',label='Query')


# Criando figura
plt.legend(loc='upper center', numpoints=1 ,ncol=3, fontsize=10, bbox_to_anchor=(0.7, 1.1))
if(len(sys.argv) >= 6):
	fig1.savefig('./results/'+str(sys.argv[6])+'.png', dpi=300)
	plt.close(fig1)
else:
	plt.show()


# Calculando delauney
# POG DIEGO - modificar isso no futuro *************************************************************************
print "Calculing delaunay"
bs.delaunay( temp , aux, sequenciasQuery)

print "| Running validation | %s"  %FileSequenciasTab
print len(sequenciasQuery)
bs.Validation(HashTabular, sequenciasQuery)

# FIM POG ******************************************************************************************************
# Fim do tempo de execucao
fim = time.time()
print "Time: ", fim-ini
