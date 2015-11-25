#!/usr/bin/python
#     Program: biosvd.py
#    Function: Implementa o metodo de svd para classificacao de proteinas
# Description: inclui: delaunay.py e sort.py - Por enquanto estamos limitado a 23 cores no plot assim so podemos ter 23 familas diferentes. Isso sera revisado(Por Thiago)
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 6


# Hierarquia do programa ******************************************
#
# ./        	Diretorio raiz
# biosvd.py 	Script principal
# README.txt	Contem instruoes de execucao
# includes  	Diretorio que armazena scripts incluidos
# tmp       	Diretorio para processamentos temporarios
# example   	Diretorio onde ficaram armazeados dados para testes
# results   	Diretorio onde deverao ser salvos os outputs	 


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


def CriarHashTab(nameFileTabFamily ):
	HashTab = {}
	Familias = []
	FileTab = open(nameFileTabFamily, 'rU')

	while True:
		row = FileTab.readline().split('\t')
		if len(row) == 1:
			break
		else:
			HashTab[ row[0] ] =  row[1]
			if row[1] not in Familias:
				Familias.append(row[1])

	FileTab.close()
	return HashTab , Familias


def Sort( sequenciasModelo, HashTabular ):

	Todasfamilias=[] #Sem repeticao
	distribuicaoFamiliasModelos = []
	SequenciasOrdenadas=[]
	FamiliaNaOrdem = [] #Ordem com que as familias aparecem no arquivo. (Com repeticao )
	seqs = []
	print "sort.py: parsing"

	for seq in sequenciasModelo:
		family = HashTabular.get(seq.id) 
		if family is not None:
			family = HashTabular.get(seq.id).rstrip() 
			FamiliaNaOrdem.append( family )

			if family not in Todasfamilias:
				Todasfamilias.append(family )

	k = -1
	for fam in Todasfamilias:
		for i in range(len(FamiliaNaOrdem)):
			if fam == FamiliaNaOrdem[i]:
				k = k+1
				seqs.append(i)
		distribuicaoFamiliasModelos.append(k)

	#gravando as sequenciasModelo na lista
	for i in seqs:
		SequenciasOrdenadas.append(sequenciasModelo[i])

	return len(Todasfamilias) ,Todasfamilias , distribuicaoFamiliasModelos , SequenciasOrdenadas





def Validation( HashTab , SeqQuery):

	filetab_prediction = open("./results/family_prediction.txt", "r")
	predictionList = []	

	tab_prediction = filetab_prediction.readline()
	tab_prediction = filetab_prediction.readline()

	while tab_prediction:
		tab_prediction = tab_prediction.split('\t')
		predictionList.append( tab_prediction[1] )

		tab_prediction = filetab_prediction.readline()

	acertos = 0
	for i in range(len(SeqQuery)):
		family = HashTab.get( SeqQuery[i].id )
		#print family,predictionList[i]
		if predictionList[i] == family:
			acertos = acertos + 1

	# Calculando o total de acertos
	porcentagem = 100 * float(acertos) / float(len(predictionList))
	print "%0.2f%% correct answers." %porcentagem
	if int(porcentagem) > 60:
		print ":)\n"
	else: 
		print ":(\n"



























os.system("clear")
try:
	mkdir("results")
except:
	print ""

# Controle do tempo de exucucao
ini = time.time()

# Declaracoes iniciais
TotalKMER = 8000;
FileNomeModelo = ''
FileNomeFamiliaModelo = ''
#FileNomeQuery = ''
FileNomeFamiliaQuery = ''


# Parametros e help
i = 1

while i < (len(sys.argv)):

	# m = modelo
	if( sys.argv[i] == '-m'):
		FileNomeModelo = sys.argv[i+1]

	# mt = modelo - arquivo tabular
	if( sys.argv[i] == '-mt'):
		FileNomeFamiliaModelo = sys.argv[i+1]

	# q = query
	if( sys.argv[i] == '-q'):
		FileQuery = sys.argv[i+1] #Contem as sequencias das query

	# qt = query - arquivo tabular
	if( sys.argv[i] == '-qt'):
		FileNomeFamiliaQuery = sys.argv[i+1]

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


# POG do THIAGO
print "Sorting Seq"
sequenciasModelo = list(SeqIO.parse( open("./example/yes.fasta", "r")  , "fasta"))
sequenciasQuery = list(SeqIO.parse( open("./example/n2.fasta", "r")  , "fasta"))

HashTabular, Familias = CriarHashTab( "Tyes.tab" )

# Numero de familas do modelo | Lista com todas as familias | Lista com a distruicao das Familias na lista | Lista com as sequencias do modelo
NumFamiliasModelo ,FamiliasModelo, DistribuicaFamiliasModelo , sequenciasModelo =  Sort(sequenciasModelo ,HashTabular )


# Clusterizacao 
print "Starting clusterization"	
# Declaracoes iniciais de clusterizacao 

print "Parsing files"
# Aqui faremos um PARSE dos arquivos .fasta para obter as sequencias a serem trabalhas
sequenciasQuery = list(SeqIO.parse(FileQuery, "fasta"))

tamanhosequenciasModelo = len(sequenciasModelo)
tamanhosequenciasQuery = len(sequenciasQuery)
NumroDeSequencias = tamanhosequenciasModelo + tamanhosequenciasQuery


aminoacidos = [ 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]
#KMERPADRAO = []
HashKMer = {}

t=0
print "Creating K-MERS"
#Gerando todos os KMERS possiveis
for i in aminoacidos:
	for j in aminoacidos:
		for k in aminoacidos:
			#KMERPADRAO.append( i+j+k )
			HashKMer[i+j+k] = t
			t = t+1

print 'Joining lists'
allSequences = []

#Jutando todas as sequencias numa so lista
for i in sequenciasModelo:
	allSequences.append(i)
for i in sequenciasQuery:
	allSequences.append(i)
	
#Apagando as lista que nao serao mais usadas
del sequenciasModelo[:]
#del sequenciasQuery[:]

#Preenchendo matriz de frequencia
mat8000x200 = bs.Kmer( allSequences  ,3 )
print "Length of the matrix: ",mat8000x200.shape

#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
Cores = ['b', 'g', 'r', 'c','m','y', 'k', 'w', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
		'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

print "SVD"
U, s, V = bs.SVD( mat8000x200 )

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


print "Creating graphic 3D"

# Salvando familias do modelo
#model = open("tmp/model.txt","w")
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
	
		x = aux[0:1,int(DistribuicaFamiliasModelo[ -1]) : tamanhosequenciasModelo ]
		y = aux[1:2,int(DistribuicaFamiliasModelo[ -1]) : tamanhosequenciasModelo ]
		z = aux[2:3,int(DistribuicaFamiliasModelo[ -1]) : tamanhosequenciasModelo ]
		F[FamiliasModelo[i]] = str(int(DistribuicaFamiliasModelo[i-1]) ) + ":" + str(tamanhosequenciasModelo)
		
	else:
		x = aux[0:1,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
		y = aux[1:2,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
		z = aux[2:3,int(DistribuicaFamiliasModelo[i-1]) :int(DistribuicaFamiliasModelo[i]) ]
		F[FamiliasModelo[i]] = str(int(DistribuicaFamiliasModelo[i-1]))+":"+str(int(DistribuicaFamiliasModelo[i]))

	ax.scatter(x, y, z, c=Cores[IDCor], marker='o', label=FamiliasModelo[IDCor])
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
plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=10, bbox_to_anchor=(0, 0))
if(len(sys.argv) >= 10):
	fig1.savefig('./results/'+str(sys.argv[10])+'.png', dpi=300)
	plt.close(fig1)
else:
	plt.show()
#model.close()

# Calculando delauney
# POG DIEGO - modificar isso no futuro *************************************************************************
print "Calculing delaunay"
bs.delaunay( temp , aux, sequenciasQuery )
#command = "python includes/delaunay.py" #envie X, Y, Z e F, alem das queries tX, tY, tZ e tF
#os.system(command)


print "| Running validation | %s"  %FileNomeFamiliaQuery
Validation(HashTabular , sequenciasQuery)
#command = "python includes/validation.py %s"  %FileNomeFamiliaQuery #envia o valor de qt
#os.system(command)
# FIM POG ******************************************************************************************************

# Fim do tempo de execucao 
fim = time.time()
print "Time: ", fim-ini
