
from random import randint, uniform, choice
import numpy as np
import math
import time
import pickle
#Anticorpo
class Anticorpo:
    def __init__(self, coordenadas = [],afinidade = 0):
        self.afinidade = afinidade
        self.coordenadas = coordenadas #paratopo
    #Gerar anticorpo aleatorio
    def gerar(self):
        #Gerar novas coordenadas para os anticorpos
        for i in range(10):    
            self.coordenadas.append(float(randint(0,1000))/1000)
    def __str__(self):
        s = 'Coordendas: {}\n Afinidade: {}\n'.format(self.coordenadas, self.afinidade)
        return s
def repertorio(numero):
    x = []
    for i in range(numero):
        anticorpo = Anticorpo([1])
        anticorpo.gerar()
        anticorpo.coordenadas = anticorpo.coordenadas[1:]
        x.append(anticorpo)
    return x
#Antigeno
class Antigeno:
    def __init__(self,coordenadas,tipo):
        self.coordenadas = coordenadas
        self.tipo = tipo #classificao (Maligno ou Benigno)
    def __str__(self):
        s = 'Coordendas: {}\n Tipo: {}\n'.format(self.coordenadas, self.tipo)
        return s
#Clonar de acordo afinidade
def clonar(pai,numero):
    filhos = []
    numero = int (math.floor(numero * pai.afinidade))
    for i in range(numero):
        filhos.append(pai)
    return filhos
#Mutacao de clones
def mutar(clones,beta,afinidade):
    coef_mut = math.exp(-beta*afinidade)
    clones2 = []
    for i in clones:
        mutado = Anticorpo([],1)
        for j in i.coordenadas:
            rand = np.random.normal(0, 0.5) #numero gaussiano
            while rand == 0:
                rand = np.random.normal(0, 0.5)
            x = j +(rand)*(coef_mut)
            mutado.coordenadas.append(x)
        clones2.append(mutado)
    return clones2
#Distancia Euclidiana -- Retorna
def distancia_euclidiana(antigeno, anticorpo):
    x = anticorpo.coordenadas[:]
    y = antigeno.coordenadas[:]
    tam = len(x)
    z = []
    for i in range(tam):
        z.append((x[i] - y[i]) ** 2)
    k = 0
    for i in range(tam):
        k += z[i]
    k = math.sqrt(k)
    return k
#Afinade normalizada -- Distancia normalizada
def normalizar(f, fmin, fmax):
    return (f - fmin)/(fmax - fmin)
#Ler arquivo
def ler(filename):
    f = open(filename, 'rU')
    total_data = []
    for line in f:
        data = line.split(',')
        data = data[:12]
        total_data.append(data)
    return total_data
#importa antigenos
def carregar(filename):
    data = ler(filename)
    antigenos = []
    for line in data:
        antigenos.append(Antigeno(line[1:12],line[0]))
    return antigenos

def main_b(printar = True):
    beta = 10
    cord = 10
    x = 500
    n_clones = 4
    n_bmcells = 5
    n_bcells = int(math.floor(n_bmcells/1.0 + 1))
    max_antigen = [0,0,0,0,0,0,0,0,0,0]
    min_antigen = range(1000,1010,1) 
    
    max_afi = -10
    min_afi = 10
    filename = "Memory.p" #Arquivo de celulas de memoria B

    arq = "wdbc.data.outcome_B_training"
    antigenos = carregar(arq) #Carregar lista de antigenos
    #colocando todas as coordenadas em reais
    for elemento in antigenos:
        for i in range(cord):
            elemento.coordenadas[i] = float(elemento.coordenadas[i])
    #Verificando maximos e minimos
    for elemento in antigenos:
        for i in range(cord):
            if elemento.coordenadas[i] > max_antigen[i]:
                max_antigen[i] = elemento.coordenadas[i]
        for i in range(cord):
            if (min_antigen[i] > elemento.coordenadas[i]):
                min_antigen[i] = elemento.coordenadas[i] 
    print "MAX: " + str(max_antigen) + "\n"+"MIN: " + str(min_antigen)
    max_a_file = open("Max_Antigen.p",'w')
    min_a_file = open("Min_Antigen.p",'w')
    pickle.dump(max_antigen,max_a_file)
    pickle.dump(min_antigen,min_a_file)
    max_a_file.close()
    min_a_file.close()
    #Normalizando todos os antigenos
    for elemento in antigenos:
        for i in range(cord):
            elemento.coordenadas[i] = normalizar(elemento.coordenadas[i],min_antigen[i],max_antigen[i])
    abr = repertorio(n_bcells)
    for i in range(x):
        at = choice(antigenos)
        max_afi = -10
        min_afi = 10
    #   print "*" * 20
    #   print at
    #   print "*" * 20
        #Calculando afinidade com antigeno
        for celula in abr:
            celula.afinidade =1/distancia_euclidiana(celula,at)
            if celula.afinidade > max_afi:
                max_afi = celula.afinidade
            if celula.afinidade < min_afi:
                min_afi = celula.afinidade
        #Normalizando afinidade
        for celula in abr:
            celula.afinidade = normalizar(celula.afinidade,min_afi,max_afi)
        #Gerando clones
        abr2 = []
        for celula in abr:
            celula2 = celula
            clone = clonar(celula,n_clones)
            #Mutar clones
            clone2 = mutar(clone,beta,celula.afinidade)
        # Afinidade dos mutados
            if len(clone2) > 1:
                for elemento in clone2:
                    elemento.afinidade = 1/distancia_euclidiana(elemento,at)
                    #print "AFINIDADE: ",elemento.afinidade
                    #print "*"*20
            
                for cell in clone2:
                    cell.afinidade = normalizar(cell.afinidade,min_afi,max_afi)
                    if cell.afinidade > celula2.afinidade:
                        celula2 = cell
                #print celula
            elif len(clone2) == 1: 
                clone2[0].afinidade = normalizar(clone2[0].afinidade,min_afi,max_afi)
                if clone2[0].afinidade >celula2.afinidade:
                    celula2 = clone2[0]
            #print celula
            abr2.append(celula2)
        abr = abr2[:]
        #Organizando celulas
        abr.sort(key=lambda cell: cell.afinidade, reverse=True)
        #Eliminandos as ruins
        abr2 = abr[:n_bmcells]
        abr  = abr2[:]
        #Gerando novas a serem testadas
        abr2 = repertorio(n_bcells - n_bmcells)
        #Juntando as antigas as novas
        abr = abr + abr2
    if printar == True:
        print "*"* 20
        print "novo conjunto"
    count = 1
    abr = abr[:n_bmcells]
    pickle_out = open(filename,"w")
    pickle.dump(abr,pickle_out)
    pickle_out.close()
    if printar == True:
        for elemento in abr:
            print count,'-',elemento
            count += 1
if __name__ == '__main__':
    START_TIME = time.time()
    main_b()
    print time.time()-START_TIME
        
        

