#!/usr/bin/python
#     Program: delauney.py
#    Function: Calcula pares de contatos usando delaunay e elimina redundancias
# Description: Recebe uma lista de posicoes x, y, z obtidas por SVD
#      Author: Diego Mariano, Thiago da Silva Correia, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 3


import numpy
from scipy.spatial import Delaunay,distance
import os


def DDelaunay( familias_modelo, matriz ):
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
		fr.write(final)

	fr.close()
	report.close()

	print "\n****************************************************** \nSuccess! \nResults at: '/results/family_prediction.txt'.\n******************************************************\n"
