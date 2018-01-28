import pickle
from ClonalgClassifier import Anticorpo, Antigeno, carregar, ler, normalizar
from ClonalgClassifier import distancia_euclidiana, main_b
from ClonalgClassifier_M import main_m
import time 

def maior_afinidade(at,test_b,test_m):
    distancia_b = 1000
    distancia_m = 1000
    for anticorpo in test_b:
        distancia = distancia_euclidiana(at,anticorpo)
        if distancia < distancia_b:
            distancia_b = distancia
    for anticorpo in test_m:
        distancia = distancia_euclidiana(at,anticorpo)
        if distancia < distancia_m:
            distancia_m = distancia
    if distancia_b < distancia_m:
        decisao = 'B'
    else:
        decisao = 'M'
    return decisao
def main():
    cord = 10

    antigenos = carregar("wdbc.data.outcome_B")
    antigenos2 = carregar("wdbc.data.outcome_M")
    filename = "Memory.p"
    filename2 = "Memory_M.p"
    pickle_read = open(filename,"r")
    pickle_read_M = open(filename2,"r")
    test_b = pickle.load(pickle_read)
    test_m = pickle.load(pickle_read_M)
    distancia_b = 1000
    distancia_m = 1000
    distancia = 0
    antigenos += antigenos2
    max_antigen = [0,0,0,0,0,0,0,0,0,0]
    min_antigen = range(1000,1010,1)
    certo = 0
    errado = 0
    total = 0
    decisao = 'B' 
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
    #Normalizando todos os antigenos
    for elemento in antigenos:
        for i in range(cord):
            elemento.coordenadas[i] = normalizar(elemento.coordenadas[i],min_antigen[i],max_antigen[i])

    
    for at in antigenos:
        decisao = maior_afinidade(at,test_b,test_m)
        if decisao == at.tipo:
            certo += 1
        else:
            errado += 1
        total += 1

    print float(certo)/float(total) * 100,"%"

        

    ''' i = 1
    for element in test_pickle:
        print str(i) +" - "
        print element
        i += 1 '''
if __name__ == "__main__":
    START_TIME = time.time()
    for i in range(5):
        main_b(False)
        main_m(False)
        main()
    print time.time()-START_TIME