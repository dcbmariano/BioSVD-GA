#     Program: sort.py
#    Function: ordena resultados
# Description: Separa por familia as sequencias usada para SVD. Eh incluso no biosvd.py
#      Author: Thiago da Silva Correia, Diego Mariano, Jose Renato Barroso, Raquel de Melo-Minardi
#     Version: 4


class Sort(object):
	"""Sort families sequences; Imput: nameFileModelo nameFileModeloFamily"""
	def __init__(self, nameFileModelo, nameFileModeloFamilia ):
		FileModelo = open(nameFileModelo, 'rU')
		FileModeloFamilia = open(nameFileModeloFamilia, 'rU')
		#FilseqModeloOrdenadas = open("tmp/SeqModeloAgrupadas.fasta", 'w')
		self.Todasfamilias=[]
		self.distribuicaoFamiliasModelos = []
		self.SequenciasOrdenadas=[]

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
				if family not in self.Todasfamilias:
					self.Todasfamilias.append(family)
			Aux = FileModeloFamilia.readline()
	
	
		#Gravando o numero total de familas do Modelo
		#FilseqModeloOrdenadas.write(str(len(self.Todasfamilias))+'\n')
		#self.distribuicaoFamiliasModelos.append(str))
		k = 0
		for fam in self.Todasfamilias:
			#FilseqModeloOrdenadas.write(str(fam)+'\n')
			for i in range(len(FamiliaNaOrdem)):
				if fam == FamiliaNaOrdem[i]:
					k = k+1
					seqs.append(i)
			#FilseqModeloOrdenadas.write(str(k)+'\n')
			self.distribuicaoFamiliasModelos.append(k)
		#print k

		#gravando as sequenciasModelo no arquivo
		for i in seqs:
			#SeqIO.write( sequenciasModelo[i], FilseqModeloOrdenadas, "fasta")
			self.SequenciasOrdenadas.append(sequenciasModelo[i])

		#del FamiliaNaOrdem[:]
		#del Todasfamilias[:]
		#del seqs[:]
		FileModeloFamilia.close()
		#FilseqModeloOrdenadas.close()

	def listas(self):
		""" Return: Number of families. List of families names. List of the families distribution; List of sequences"""
		return len(self.Todasfamilias) ,self.Todasfamilias , self.distribuicaoFamiliasModelos ,self.SequenciasOrdenadas
