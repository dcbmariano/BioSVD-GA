#!/usr/bin/python
#     Program: validation.py
#    Function: calcula a porcentagem de acertos
# Description: eh incluido por biosvd.py
#      Author: Diego Mariano, Thiago da Silva Correia, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 2

# Imports
import sys

# Recebendo argumentos
#tab_query_file = sys.argv[1]
#tab_query_file = "../example/amostra.tab" # para rodar diretamente do arquivo, descomente esta linha

class Validation(object):

	def __init__(self, tab_query_file ):

		tab_query = open(tab_query_file).readlines()
		tab_prediction = open("results/family_prediction.txt").readlines()
		#tab_prediction = open("../results/family_prediction.txt").readlines() # rodar diretamente do arquivo, descomente esta linha

		p = []
		q = []

		# Log de erros
		erro = open("tmp/log_erros.txt","w")

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
