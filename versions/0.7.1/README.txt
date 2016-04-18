# BIOSVD 
# VERSION: 0.7.1 ALPHA

# Documentation
# Language: PT-BR

# Hierarquia do programa ******************************************
#
# ./        	Diretorio raiz
# biosvd.py 	Script principal
# README.txt	Contem instruoes de execucao
# includes  	Diretorio que armazena scripts incluidos
# tmp       	Diretorio para processamentos temporarios
# example   	Diretorio onde ficaram armazeados dados para testes
# results   	Diretorio onde deverao ser salvos os outputs	


# Parametros ******************************************
- q             Nome do arquivo query (formato fasta)
- qt            Arquivo tabular para query (formato tab)
- m             Nome do arquivo modelo (formato fasta)
- mt            Arquivo tabular para query (formato tab)
- g             Salvar grafico (insira o nome do arquivo sem extensao)
- p				Nome do arquivo para salvar o plot do posto da matriz(insira o nome do arquivo sem extensao)



# Uso simples ******************************************

Para rodar o script sem a vaildacao cruzada
python example.py -A -m example/yes.fasta -mt example/YES.tab -q example/no.fasta -qt example/NO.tab -g plotClus -p plotPosto

Para rodar o script com a validaao cruzada
Opcao de entrada 1
python Cross.py -A -k 10 -m example/yes.fasta -mt example/YES.tab -q example/no.fasta -qt example/NO.tab -g plotClus -p plotPosto

Opcao de entrada 2
python Cross.py -B -g plotClus -p plotPosto -k Kfold [arquivos fasta ] 
