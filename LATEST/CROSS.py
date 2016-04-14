#!/usr/bin/python
#     Program: experimento.py
#    Function: Implementa o experimento de Dobson com o metodo biosvd
# Description: 
#      Author: Diego Mariano, Thiago da Silva Correia, Jose Renato Barroso, Raquel Cardoso de Melo-Minardi
#     Version: 1


import sys
from BioSVD import BioSVD 
from Bio import SeqIO
from numpy import * 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from random import randint
from scipy.spatial import Delaunay,distance


# Parametros e help
i = 1
familias = []
seqs = ''
inicio = 0
fim = 0

# Percorre os parametros enviados ########################################################
while i < (len(sys.argv)):


	# Help ###############################################################################
	if sys.argv[i] == '-h' or sys.argv[i] == '--help':
		print "\n*** BIOSVD ***\nSyntax: \npython experimento.py [arquivos.fasta]\n"
		sys.exit()

	# End Help ###########################################################################


	# Cada arquivo fasta contem uma familia diferente, logo so precisa do fasta pra rodar
	arquivo = sys.argv[i]
	familia_completo = arquivo.split("/") # caso envie a partir de varias pastas
	familia_split = familia_completo[len(familia_completo)-1].split(".") # remove a extensao
	familia = familia_split[0]

	# Lendo CADA sequencia com biopython
	seq = list(SeqIO.parse( open(arquivo, "r")  , "fasta"))
	if seqs != '':
		seqs = seqs + seq
	else:
		seqs = seq

	total_seqs = len(seq)

	fim = inicio + total_seqs 

	# o que familias armazena - Ex.:
	# hidrolase, 0, 100
	# Define uma cor aleatoria 
	numero_aleatorio = randint(100000,999999)
	familias.append([familia,inicio,fim,numero_aleatorio])

	# incrementa o contador
	i += 1
	inicio = fim


# Criando matriz de k-mers ##################################################################
matkmer = biosvd.Kmer(seqs,3)

# SVD #######################################################################################
print "SVD"
U, s, V = linalg.svd(matkmer)
SK =  diag(s) # Normlizando a matriz S
Temp = SK/sum(s) 
SN = diag( Temp )

# SVD | Visualizando posto
fig2 = plt.figure()
plt.axis([0, 20, 0, 0.5]) #Aqui limito os eixos X e Y do plot
plt.plot(SN.T)
fig2.savefig('posto.png', dpi=300)
plt.close(fig2)

# SVD | Criando matriz de posicionamento
posto = 3 # usa posto 3 por default
SK = SK[0:posto,0:posto]
Vaux = V.transpose()
VK = Vaux[:,:posto]
aux = dot(SK , VK.transpose() )

# SVD | Visualizando matriz de posicionamento
fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d') # Requer from mpl_toolkits.mplot3d import Axes3D

# Preenche as posicoes x, y e z
for i in range(fim):
	x = aux[0:1, i:i+1]
	y = aux[1:2, i:i+1]
	z = aux[2:3, i:i+1]

	# Definindo cores
	for j in familias:
		if i >= j[1] and i < j[2]:
			cor = "#%d" %(j[3])

	ax.scatter(x, y, z, c=cor, marker='o', s=100, label=j[0]) # plota o ponto

plt.show()
plt.close(fig1)

# Invertendo a matriz aux -> pontos (transpondo)
pontos = []
tam = len(aux[0])
for i in range(tam):
	pontos.append([aux[0][i],aux[1][i],aux[2][i]])



# Validacao cruzada ######################################################################

# seleciona 5 ou 10 % -> regra: quantidade proporcionalmente igual para todos
# O teste sera realizado com 5 ou 10% da base, enquanto o treino eh feito com 90 ou 95%
vc = 10; #5 ou 10 (10 default) 

# CRIA OS GRUPOS
# Atencao: ele distribui igualmente os elementos em grupos 
# a quantidade de grupos eh definida na variavel vc
# Ele le cada elemento (variavel tam) e armazena no grupo 0
# o proximo elemento ele armazena no grupo 1, o proximo no 2
# e assim sucessivamente
# os grupos estao armazenados em fold e as familias em fold_familia 
j = 0
fold = []
fold_familia = []

# Crie os grupos em branco (nao consegui resolver de outra maneira :[ )
for i in range(vc):
	fold.append([])
	fold_familia.append([])

# Percorre cada elemento
# tam armazena tamanho do array pontos
for i in range(tam):

	# quando chegar ao ultimo grupo (=vc), retorne o valor a zero
	if j == vc:
		j = 0

	# Cada elemento sera inserido no grupo j
	fold[j].append(pontos[i])

	for f in familias:
		if i >= f[1] and i < f[2]: # f[0] armazena o nome da familia, f[1] inicio em tam, f[2] fim em tam
			fold_familia[j].append(f[0]) # esse codigo merece ganhar o nobel da programacao

	j += 1


# Arquivo de resultado
w = open("result.txt","w")

acertos_teste = 0
acertos_erros_teste = 0
acertos_total = 0
acertos_erros_total = 0

# Agora que os grupos estao separados vamos CLASSIFICAR
# O teste sera repetivo 5 ou 10 vezes
for i in range(vc):

	print "\nREALIZANDO TESTE %d\n" %(i)

	# Gravando resultados
	titulo = "::Teste %d\n" %(i)
	menu = "Proteina_Grupo\tAcerto\tPredicao\tOriginal\tVotos\n"
	w.write(titulo)
	w.write(menu)

	treino = fold[:] #faz a copia completa de fold :) -> se eu fizer apenas treino = fold ele copia apenas a posicao na memoria
	teste = treino.pop(i) #elimina i de treino e retorna para teste

	treino_familia = fold_familia[:]
	teste_familia = treino_familia.pop(i)

	# Para cada elemento proteina no grupo i
	proteina_id = 0 # necessario pois proteina eh um objeto e nao um inteiro
	for proteina in fold[i]:

		print "Analisando proteina %d do grupo %d." %(proteina_id,i)
		
		
		# roda delaunay
		# Aqui faco o delauney de um dos elementos da base teste + toda base de treino
		teste_atual = [] #zero o valor de teste atual
		teste_atual.append(proteina) # insiro o valor testado primeiro

		# isso eh necessario pois ha listas dentro de listas, tenho que colocar todo mundo no mesmo nivel
		for cada_grupo in treino: 
			for elemento in cada_grupo:
				teste_atual.append(elemento)

		# O mesmo deve ser feito em familia -> assim familia_teste_atual equivale a p[1] ou p[0]
		familia_teste_atual = []
		familia_teste_atual.append(fold_familia[i][proteina_id]) # elemento que devera ser predito :) 

		for cada_grupo in treino_familia: 
			for elemento in cada_grupo:
				familia_teste_atual.append(elemento) 
		

		# ----------------------------------- CLIMAX ------------------------------------- #
		delaunay = Delaunay(teste_atual) 

		# Extraindo os pares
		pares = []
		for par in delaunay.vertices: 
			for k in range(4):
				for j in range(k+1,4):
					pares.append([par[k],par[j]])
		#print "Total of pairs: ",len(pares)

		# Remove duplicacoes 
		pares_unicos = []
		for x in pares:
			if x not in pares_unicos:
				pares_unicos.append(x)
		#print "Total of unique pairs: ",len(pares_unicos)

		# Busca contatos
		votos = {}
		for f in familias:
			votos[f[0]] = 0 # zera o contador de votos

		for p in pares_unicos:
			if p[0] == 0:
				print "Proteina %d (familia %s) faz contato com %d (familia %s)." %(proteina_id,familia_teste_atual[0],p[1],familia_teste_atual[p[1]])
				
				# ALGORITMO DE VOTACAO - usado pra predizer a familia
				#para cada familia
				for f in familias:
					# nome da familia eh armazenado em f[0]
					if f[0] == familia_teste_atual[p[1]]:
						votos[f[0]] += 1


			if p[1] == 0:
				print "Proteina %d (familia %s) faz contato com %d (familia %s)." %(proteina_id,familia_teste_atual[0],p[0],familia_teste_atual[p[0]])

				# ALGORITMO DE VOTACAO - usado pra predizer a familia
				#para cada familia
				for f in familias:
					# nome da familia eh armazenado em f[0]
					if f[0] == familia_teste_atual[p[0]]:
						votos[f[0]] += 1
		
		print votos

		# Gravando resultados
		proteina_grupo = "%d_%d" %(proteina_id,i)
		acerto = 0
		predicao = max(votos,key=votos.get)
		original = familia_teste_atual[0]
		if predicao == original:
			acerto = 1

		acertos_teste += acerto
		acertos_erros_teste += 1
		acertos_total += acerto
		acertos_erros_total += 1

		linha = "%s\t%d\t%s\t%s\n" %(proteina_grupo, acerto, predicao, original)
		print linha
		w.write(linha)

		proteina_id += 1

	# Calculando acuracia
	acuracia_teste = float(acertos_teste)/float(acertos_erros_teste)
	acertos_teste = 0 # zera variaveis
	acertos_erros_teste = 0 # zera variaveis

	print_acuracia = "\nAcuracia do teste: %0.2f\n\n" %(acuracia_teste)
	w.write(print_acuracia)
	print print_acuracia
	w.write("--------------------------------------------------------------\n\n")

acuracia_total = float(acertos_total)/float(acertos_erros_total)
print "acertos total: "
print acertos_total

print "acertos erros: "
print acertos_erros_total

print "acuracia_total: "
print acuracia_total

print_acuracia_total = "\nAcuracia media: %0.2f\n\n" %(acuracia_total)
print print_acuracia_total
w.write(print_acuracia_total)
print "\n--------------------------------------------------------------\n"
print "Obrigado por usar BioSVD.\nFim da execucao.\n"
print "--------------------------------------------------------------\n"
w.close()
