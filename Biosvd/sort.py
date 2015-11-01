#     Program: sort.py
#    Function: ordena resultados
# Description: Separa por familia as sequencias usada para SVD. Eh incluso no biosvd.py
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel de Melo-Minardi
#     Version: 4


from Bio import SeqIO


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


	#gravando as sequenciasModelo na lista
	for i in seqs:
		SequenciasOrdenadas.append(sequenciasModelo[i])

	FileModeloFamilia.close()

	return len(Todasfamilias) ,Todasfamilias , distribuicaoFamiliasModelos , SequenciasOrdenadas
