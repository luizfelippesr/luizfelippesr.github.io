import random

with open('texto_secreto.txt','r') as f:
    texto_secreto = f.read()

linhas = texto_secreto.splitlines()
palavras = []
x = []
for l in linhas:
    x.append(len(l.split(' ')))
    palavras += l.split(' ')

segredo = 'xyi'

for i, linha_secreta in enumerate(linhas):
    print('{0:4d}'.format(i),sep='\b',flush=True)
    # Define o arquivo que contém uma linha secreta
    f = open('secreto/arquivo{0}.txt'.format(i+1),'w+')

    # Define a localização da linha secreta no arquivo
    isecreta = random.randint(0,len(linhas))
    for ilinha in range(len(linhas)):
        if isecreta == ilinha:
            # Se for a linha sorteada, imprime a linha secreta
            f.write(linha_secreta + segredo+'\n')
        else:
            # Caso contrário, imprime entre 9 e 11 palavras
            n = random.randint(9,11)
            for j in range(n):
                f.write(random.choice(palavras)+' ')
            f.write('\n')
    f.close()
