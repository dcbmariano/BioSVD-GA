#!/usr/bin/env python
#
#     Program: biosvd.py
#    Function: Implementa o metodo de svd para classificacao de proteinas
# Description: inclui: delaunay.py e sort.py - Por enquanto estamos limitado a 23 cores no plot assim so podemos ter 23 familas diferentes. Isso sera revisado(Por Thiago)
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 7.01


from numpy import linalg as LA
from scipy.sparse.linalg import svds
from Bio import SeqIO
import itertools
import numpy as np
from scipy.spatial import Delaunay,distance
import multiprocessing as mp
import matplotlib.pyplot as plt
import math
import os
from functools import partial
from random import *
import sys

def SVD (  matriz , K , namePlotPosto):

	U, s, V = svds( matriz )
	V = V.transpose()
	tam = len(s)
	s[:] =  s[::-1] #Ordenando de forma decrescente
	S = np.diag(s)	#Criando a matrix S cuja colunas sao formandas pelo vetor s
	SK = S[0:K,0:K]
	VK = V[:,:K]
	aux = np.dot(SK , VK.transpose() )
	Posto( s,namePlotPosto )
	return aux, U
 
def Kmer( sequencias, k ):
	print "Creating K-MERS"
	aminoacidos = [ 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]
	HashKMer = {}
	t = 0
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

#Para exibir o posto para passar o namePlotPosto vazio
def Posto( vetorS , namePlotPosto ):
	S = vetorS*(1/np.sum(vetorS))	#Dividindo pela soma do vetor s
	fig2 = plt.figure()
	fig2.suptitle('Posto')
	#Aqui limito os eixos X e Y do plot
	#plt.axis([0, tam , 0, 0.5])
	plt.plot(S)
	if(len(namePlotPosto) >0):
		fig2.savefig('results/'+namePlotPosto+'.png', dpi=300)
		plt.close(fig2)
	else:
		plt.show()


def delaunay( familias_modelo, matriz , sequenciasQuery, HastabularQuey):
	# definindo familias
	opcao_entrada = '-A'
	f = {}
	for i in range(len(familias_modelo)):
		if i % 2 == 0:
			aux = familias_modelo[i].strip()
		if i % 2 == 1:
			f[aux] = familias_modelo[i].strip()
			ultimo = familias_modelo[i].split(":")
			ultimo_elemento_modelo = int(ultimo[1])

	# Recebendo todos os dados (variaveis x, y e z)
	x = matriz[0]
	y = matriz[1]
	z = matriz[2]

	# Resultado final
	fr = open("./results/family_prediction.txt","w")
	#fID = open("./results/ID_prediction.txt","w")
	fr.write("protein_id (query)\tfamily predicted\tfamily correct\n")


	# Predicao aqui
	print "Determining contacts"

	# Relatorio de percentual
	report = open("./results/report_contacts.txt","w")

	# Modificacao aqui
	# delauney deve ser executado com o modelo e cada proteina individualmente
	coord = []

	# Unindo matrizes x y z => APENAS MODELO
	for i in range(len(x)):
		coord.append([float(x[i]),float(y[i]),float(z[i])])

	# IMPORTANTE
	# Inseriremos um elemento por vez na lista coord e calculamos o delauney - matriz x eh usada como exemplo
	jend = len(x)

	for i in range(ultimo_elemento_modelo,jend):

		# UM POR VEZ  
		lista_atual = coord[:ultimo_elemento_modelo] + [coord[i]]
		tam_lista_atual = len(lista_atual)

		# Calculando os tetraedros ********************************************* PONTO MAIS IMPORTANTE AQUI =>
		delaunay = Delaunay(lista_atual)

			# Extraindo os pares ********************************************
		pares = []
		for par in delaunay.vertices: 
			for k in range(4):
				for j in range(k+1,4):
					pares.append([par[k],par[j]])

		  #print "Total of pairs: ",len(pares)

		  # Remove duplicacoes ********************************************
		pares_unicos = []
		for x in pares:
			if x not in pares_unicos:
				pares_unicos.append(x)

	# Validando o teste - proteinas testes estao armazenadas em pares_unicos[da ultima posicao do modelo ate a ultima posicao da lista coord]

		# Contador de elementos por familia (cef)
		# Deve ser zerado para cada elemento testado
		cef = {}
		for k in f:
			cef[k] = 0

		# Analisamos cada par 
		for j in pares_unicos:
			# Verificamos se o elemento esta na posicao 0
			if j[0] == tam_lista_atual-1:
				# para cada familia F armazenada na lista f
				for F in f:
					posicao = f[F].split(":")
					if j[1] >= int(posicao[0]) and j[1] <= int(posicao[1]):
						cef[F] += 1

			# AGORA PRECISAMOS CONFERIR SE O VALOR TESTADO ESTA NA SEGUNDA POSICAO DO PAR
			# Verificamos se o elemento esta na posicao 1
			if j[1] == tam_lista_atual-1:
			# para cada familia F armazenada na lista f
				for F in f:
					posicao = f[F].split(":")
					if j[0] >= int(posicao[0]) and j[0] <= int(posicao[1]):
						cef[F] += 1

		# determina a familia com base em quem tem o maior numero de contatos
		familia_atual = max(cef, key=cef.get)

		# Relatorio de percentual
		report_txt = ">position %d [%d:%d]\n" %(i,ultimo_elemento_modelo,len(coord))
		report.write(report_txt)
		report.write(str(cef))
		report.write("\n")
		  
		# calcula o id da protein removendo a quantidade do modelo
		protein_id = i - ultimo_elemento_modelo + 1
		if opcao_entrada ==  '-A':
			ID = sequenciasQuery[protein_id-1].id.split("|")[2]
		else:
			ID = sequenciasQuery[protein_id-1].id
		family_correct = HastabularQuey.get(ID)
		final = "%d\t%s\t%s\n" %(protein_id, familia_atual, str(family_correct))
		#temp = 	"%s\t%s\n" %(str(sequenciasQuery[protein_id-1].id),familia_atual)
		#fID.write(temp)
		fr.write(final)

	fr.close()
	report.close()

	print "\n****************************************************** \nSuccess! \nResults at: '/results/family_prediction.txt'.\n******************************************************\n"




def Validation( HashTab , SeqQuery ):
	filetab_prediction = open("./results/family_prediction.txt", "r")
	predictionList = []	
	opcao_entrada = '-A'

	tab_prediction = filetab_prediction.readline()
	tab_prediction = filetab_prediction.readline()

	while tab_prediction:
		tab_prediction = tab_prediction.split('\t')
		predictionList.append( tab_prediction[1].rstrip() )

		tab_prediction = filetab_prediction.readline()

	acertos = 0
	proteinas_sem_familias = 0
	erros = 0
	for i in range(len(SeqQuery)):
		if opcao_entrada ==  '-A':
			ID = SeqQuery[i].id.split("|")[2]
		else:
			ID = SeqQuery[i].id
		family = HashTab.get(ID)
		if family is not None:
			family = family.rstrip()
		
			if i <len(predictionList) and predictionList[i] == family:
				acertos = acertos + 1
			else:
				erros = erros +1
		else:
			proteinas_sem_familias = proteinas_sem_familias + 1

	# Calculando o total de acertos
	porcentagem = 100*float(acertos) / float(len(predictionList) - proteinas_sem_familias )
	print "%0.2f%% correct answers." %porcentagem
	if int(porcentagem) > 60:
		print ":)\n"
	else: 
		print ":(\n"
	return acertos,erros



#Agrupa as sequencias por familia
def Sort( sequenciasModelo, HashTabular ):
	opcao_entrada= '-A'
	Todasfamilias=[] #Sem repeticao
	distribuicaoFamiliasModelos = []
	SequenciasOrdenadas=[]
	FamiliaNaOrdem = [] #Ordem com que as familias aparecem no arquivo. (Com repeticao )
	seqs = []

	for seq in sequenciasModelo:
		if opcao_entrada =='-A':
			ID = seq.id.split("|")[2]
		else:
			ID = seq.id
		family = HashTabular.get(ID)
		if family is not None:
			family = family.rstrip() 
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



#Criar uma tabela Hash onde os indices e o ID das sequencias
def CriarHashTab( nameFileTabFamily ):
	HashTab = {}
	Familias = []
	i = 0
	FileTab = open(nameFileTabFamily, 'rU')
	cont =0
	while True:
		cont = cont +1
		row = FileTab.readline().split('\t')
		if len(row) == 1:
			break
		else:
			HashTab[ row[0].rstrip() ] =  row[1].rstrip()
			if row[1] not in Familias:
				Familias.append(row[1].rstrip())

	FileTab.close()
	return HashTab , Familias

def Cosseno( A, B ):
	dot = A* B 
	normA = LA.norm(A)
	normB = LA.norm(B)

	cos = float( dot/( normA * normB) )
	return cos

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


def CrossValidation( sequencias, hashTab , kfold ):

	'''
	temp = []
	while len(temp) < 1000:
		try:
			i =randint(0, len(sequencias))
			temp.append(sequencias[i])
			del sequencias[i]
		except:
			print ""

	sequencias = temp
	print "TAMANHO SEQ: ", len(sequencias)
	'''
	numeroSequencias = len(sequencias)
	numero_Seq_por_grupos= numeroSequencias/kfold
	ruido = 0.0
	

	#SlateBlue, MediumVioletRed, DarkOrchid,DeepSkyBlue,DarkRed,OrangeRed,Teal,
	#Lime,DarkGoldenrod,PaleTurquoise,Plum,LightCoral,CadetBlue,DarkSeaGreen,PaleGoldenrod,RosyBrown
	Cores = ['b', 'g', 'r', 'c','m','y', 'k', 'w', '#6A5ACD', '#C71585','#9932CC','#8B0000','#FF4500',
	'#008B8B','#00FF00','#B8860B','#E0FFFF','#DDA0DD' ,'#F08080' ,'#5F9EA0','#8FBC8F','#EEE8AA','#BC8F8F']

	#Caso o numero de sequencias nao seja multiplo do numero numero_Seq_por_grupos * kfold
	#As ultima sequencias serao ignoradas

	print "Creating K-fold"
	for j in range( kfold ):
		print j
		seq = []
		while len(seq) < numero_Seq_por_grupos :
			try:
				i =randint(0, len(sequencias))
				seq.append(sequencias[i])
				del sequencias[i]
			except:
				print ""
		SeqIO.write(seq, open("./results/grupo"+str(j)+".fasta", "w"), "fasta")
		del seq[:]
	'''
	print "number sequence: %d \nk-fold: %d\nsequences by groups: %d" %( len(sequencias), kfold, numero_Seq_por_grupos)

	# Aqui vamos seperar as sequencias modelo e as Query
	for i in range(numero_Seq_por_grupos ):
		nomePlotClusterizacao = "Clust Group: " + str(i)
		print "Grupo: %d" %int(i+1)

		print "Sorting Model Seq"
		sequenciasModel = sequencias[-numero_Seq_por_grupos:]
		NumFamiliasModelo,FamiliasModelo,DistribuicaFamiliasModelo,sequenciasModel=Sort( sequenciasModel ,hashTab)

		#Preenchendo matriz de frequencia
		matFrequencia = Kmer( sequencias  ,3 )
		print "Length of the matrix: ",matFrequencia.shape
		print "SVD"
		aux ,T = SVD( matFrequencia, 3 , "Posto Grupos "+str(i) )

		numeroSeqModel = len(sequenciasModel)

		print "Creating graphic 3D"

		fig1 = plt.figure()
		ax = fig1.add_subplot(111, projection='3d')

		F = {}
		tF = {}
		IDCor = 0
		temp = []


		#print "Model"
		for l in range(len(DistribuicaFamiliasModelo)):
			if l == 0:#Primeira Familia Modelo
				x = aux[0:1,0:int(DistribuicaFamiliasModelo[0]) ]
				y = aux[1:2,0:int(DistribuicaFamiliasModelo[0]) ]
				z = aux[2:3,0:int(DistribuicaFamiliasModelo[0]) ]
				F[FamiliasModelo[l]] = "0:"+str(int(DistribuicaFamiliasModelo[0]) )
		
			elif(l == len(DistribuicaFamiliasModelo)-1 ):#Ultima familia Modelo
	
				x = aux[0:1,int(DistribuicaFamiliasModelo[ -1]) : numeroSeqModel ]
				y = aux[1:2,int(DistribuicaFamiliasModelo[ -1]) : numeroSeqModel ]
				z = aux[2:3,int(DistribuicaFamiliasModelo[ -1]) : numeroSeqModel ]
				F[FamiliasModelo[l]] = str(int(DistribuicaFamiliasModelo[l-1]) ) + ":" + str(numeroSeqModel)
		
			else:
				x = aux[0:1,int(DistribuicaFamiliasModelo[l-1]) :int(DistribuicaFamiliasModelo[l]) ]
				y = aux[1:2,int(DistribuicaFamiliasModelo[l-1]) :int(DistribuicaFamiliasModelo[l]) ]
				z = aux[2:3,int(DistribuicaFamiliasModelo[l-1]) :int(DistribuicaFamiliasModelo[l]) ]
				F[FamiliasModelo[l]] = str(int(DistribuicaFamiliasModelo[l-1]))+":"+str(int(DistribuicaFamiliasModelo[l]))

			ax.scatter(x, y, z, c=Cores[IDCor], marker='o', label=FamiliasModelo[IDCor], s = 100 )
			IDCor = IDCor +1
			temp.append(FamiliasModelo[l])
			temp.append(str(F[FamiliasModelo[l]]))

	
		#print "Query"
		tx = aux[0:1, numeroSeqModel +1: ]
		ty = aux[1:2, numeroSeqModel +1: ]
		tz = aux[2:3, numeroSeqModel +1: ]
		tF['Query family'] = str(numeroSeqModel+1)+":"+str(numeroSequencias )
		ax.scatter(tx, ty, tz, c='#000000', marker='*',label='Query', s = 100 )

		# Criando figura
		plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=10, bbox_to_anchor=(0, 0))
		if(len(nomePlotClusterizacao) > 0):
			fig1.savefig('results/'+nomePlotClusterizacao+'Grupo'+str(i) +'.png', dpi=300)
			plt.close(fig1)
		else:
			plt.show()

		# Calculando delauney
		sequenciasQuery = sequencias[0:-numero_Seq_por_grupos]
		print "ABBBBBBBBBBBBBB:",len(sequenciasQuery)
		#print "Calculing delaunay"
		delaunay( temp, aux , sequenciasQuery, hashTab )

		#print "| Running validation |" 
		ruido = (1 - Validation( hashTab, sequenciasQuery) ) + ruido

		seqTemp = sequencias[ 0:numero_Seq_por_grupos ]
		del sequencias[ 0:numero_Seq_por_grupos ]
		sequencias = sequencias + seqTemp

	print "AC: " ,float(ruido/numeroSequencias)


	return'''










