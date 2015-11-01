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
'''
def Sort( nameFileModelo, nameFileModeloFamilia ):
	FileModelo = open(nameFileModelo, 'rU')
	FileModeloFamilia = open(nameFileModeloFamilia, 'rU')
	#FilseqModeloOrdenadas = open("tmp/SeqModeloAgrupadas.fasta", 'w')
	Todasfamilias=[]
	distribuicaoFamiliasModelos = []
	SequenciasOrdenadas=[]

	print "sort.py: parsing"
	# Aqui faremos um PARSE dos arquivos .fasta para obter as sequencias a serem trabalhas
	sequenciasModelo = list(SeqIO.parse(FileModelo, "fasta"))
	FileModelo.close()
	seqs = []

	FamiliaNaOrdem = [] #Ordem com que as familias aparecem no arquivo. (Com repeticao )

	#Modelo
	Aux = FileModeloFamilia.readline()
	Aux = FileModeloFamilia.readline()
	while ( Aux ):
		temp = Aux.split('\t')
		tempA = temp[7].split(',')
		i = tempA[0].rstrip().find("Glycosyl hydrolase")
		# Faco a busca apartir da posicao em que se encotra o Glycosyl hydrolase
		j = tempA[0].rstrip().find("family",i, len(tempA[0].rstrip()))
		family = tempA[0].rstrip()[i:j-1]
		if(len(family)>0):
			FamiliaNaOrdem.append( family )
			if family not in Todasfamilias:
				Todasfamilias.append(family)
		Aux = FileModeloFamilia.readline()

	#Gravando o numero total de familas do Modelo
	k = 0
	for fam in Todasfamilias:
		#FilseqModeloOrdenadas.write(str(fam)+'\n')
		for i in range(len(FamiliaNaOrdem)):
			if fam == FamiliaNaOrdem[i]:
				k = k+1
				seqs.append(i)
		distribuicaoFamiliasModelos.append(k)


	#gravando as sequenciasModelo no arquivo
	for i in seqs:
		SequenciasOrdenadas.append(sequenciasModelo[i])

	FileModeloFamilia.close()
	#FilseqModeloOrdenadas.close()

	""" Return: Number of families. List of families names. List of the families distribution; List of sequences"""
	return len(Todasfamilias) ,Todasfamilias , distribuicaoFamiliasModelos ,SequenciasOrdenadas


def Validation( tab_query_file ):

	tab_query = open(tab_query_file).readlines()
	tab_prediction = open("family_prediction.txt").readlines()
		
	p = []
	q = []

	# Log de erros
	erro = open("log_erros.txt","w")

	# Preenchendo a matriz da query
	for i in range(len(tab_query)):

		# Ignore a primeira linha: cabecalho
		if i != 0:

			# Quebra com base na tabulacao
			linha = tab_query[i].split("\t")

			# Quebra por espacamento para retirar apenas o numero (esta na posicao 7)
			try:
				num = linha[7].split(" ")
			except:
				fail_msg = "Fail to split %d. Validation.py will Ignore and continue.\n" %i
				fail_msg2 = num
				erro.write(fail_msg)
				erro.write(str(fail_msg2))

			# Insere na lista q 
			try:
				q.append(num[2].strip())
			except:
				fail_msg = "Fail to split %d. Validation.py will Ignore and continue.\n" %i
				fail_msg2 = num
				erro.write(fail_msg)
				erro.write(str(fail_msg2) )



	# Preenchendo a matriz da prediction
	for i in range(len(tab_prediction)):

		# Ignore a primeira linha: cabecalho
		if i != 0:

			# Quebra com base na tabulacao
			linha = tab_prediction[i].split("\t")

			# Quebra por espacamento para retirar apenas o numero (esta na posicao 1)
			try:
				num = linha[1].split(" ")
			except:
				fail_msg =  "Fail to split %d. Validation.py will Ignore and continue.\n" %i
				fail_msg2 =  num
				erro.write(fail_msg)
				erro.write(str(fail_msg2))

			# Insere na lista p
			try:	
				p.append(num[2].strip())
			except:
				fail_msg =  "Fail to split %d. Validation.py will Ignore and continue.\n" %i
				fail_msg2 = num
				erro.write(fail_msg)
				erro.write(str(fail_msg2))


	# Comparando
	total_elementos = len(p)
	acertos = 0

	for i in range(total_elementos):
		try:
			#print q[i], p[i]
			if q[i] == p[i]:
				acertos += 1
		except:
			fail_msg =  "Fail to Compare %d. Validation.py will Ignore and continue.\n" %i
			erro.write(str(fail_msg2))

	# Calculando o total de acertos
	porcentagem = 100 * acertos / total_elementos

	print "%i%% correct answers." %porcentagem
	if int(porcentagem) > 60:
		print ":)\n"
	else: 
		print ":(\n"
	erro.close()

def AplciarDelaunay ( familias_modelo , matriz ):
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
	# Exemplo input
	# x = [-4.26841,-7.09907,-6.85084,-7.05957,-7.78660,-6.52208,-7.10390,-7.06432,-5.21143,-8.08927]
	# y = [-2.84498,-6.00829,-6.12811,-5.66237,-6.58640,-5.44261,-5.23417,-5.16905,-4.44820,-6.40365]
	# z = [0.81944,1.33445,1.00721,1.30053,1.94627,1.01903,0.70142,0.93283,1.69986,1.87179]
	x = matriz[0]
	y = matriz[1]
	z = matriz[2]

	# Resultado final
	fr = open("family_prediction.txt","w")
	fr.write("protein_id (query)\tfamily predicted\n")


	# Predicao aqui
	print "Determining contacts"

	# Relatorio de percentual
	report = open("report_contacts.txt","w")

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
		fr.write(final)

	fr.close()
	report.close()

	print "\n****************************************************** \nSuccess! \nResults at: 'results/family_prediction.txt'.\n******************************************************\n"
'''
