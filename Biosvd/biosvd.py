#!/usr/bin/env python
# -*- coding: latin-1 -*-


from numpy import linalg as LA
from Bio import SeqIO
import itertools
import numpy as np
from scipy.spatial import Delaunay,distance
import math
import os

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

def delaunay( familias_modelo, matriz , sequenciasQuery ):
	# definindo familias
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
	fID = open("./results/ID_prediction.txt","w")
	fr.write("protein_id (query)\tfamily predicted\n")


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

		  #print "Total of unique pairs: ",len(pares_unicos)

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

		final = "%d\t%s\n" %(protein_id,familia_atual)
		temp = 	"%s\t%s\n" %(str(sequenciasQuery[protein_id-1].id),familia_atual)
		fID.write(temp)
		fr.write(final)

	fr.close()
	report.close()

	print "\n****************************************************** \nSuccess! \nResults at: '/results/family_prediction.txt'.\n******************************************************\n"




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



#Agrupa as sequencias por familia
def Sort( sequenciasModelo, HashTabular ):

	Todasfamilias=[] #Sem repeticao
	distribuicaoFamiliasModelos = []
	SequenciasOrdenadas=[]
	FamiliaNaOrdem = [] #Ordem com que as familias aparecem no arquivo. (Com repeticao )
	seqs = []
	print "sort.py: parsing"

	for seq in sequenciasModelo:
		family = HashTabular.get(seq.id).rstrip() 
		if family is not None:
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



#Criar uma tabela Hash onde os indices é o ID das sequencias
def CriarHashTab( nameFileTabFamily ):
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

#Calculates euclidean distances
def dist_euclidiana(p1,p2):
	return np.linalg.norm(p1-p2)


def Euclidiana( NumroDeSequencias ):
	matProteinasProximas = np.zeros(shape=( NumroDeSequencias ,5 ) )
	distancias = np.zeros(shape=( NumroDeSequencias , NumroDeSequencias ) )

	
	MatPosicao = np.loadtxt( 'MatrizPosicoes.txt', delimiter=';' )
	dimensao = MatPosicao.shape[0]


	for i in range(NumroDeSequencias):
		for j in range(NumroDeSequencias):
			if  j > i:
				p1 = MatPosicao[:dimensao,i:i+1]
				p2 = MatPosicao[:dimensao,j:j+1]
				distancias[i][j] = dist_euclidiana(p1,p2)











