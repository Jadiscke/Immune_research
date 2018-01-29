import pickle
from ClonalgClassifier import Anticorpo, Antigeno, carregar, ler, normalizar
from ClonalgClassifier import distancia_euclidiana, main_b
from ClonalgClassifier_M import main_m
import time 


def maior_afinidade(at,test_b,test_m):
    #retorna o chute da antigeno se baseando na celula de memoria mais proxima
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
    percentuais = []
    antigenos = carregar("wdbc.data.outcome_B")
    antigenos2 = carregar("wdbc.data.outcome_M")
    filename = "Memory.p"
    filename2 = "Memory_M.p"
    pickle_read = open(filename,"r")
    pickle_read_M = open(filename2,"r")
    test_b = pickle.load(pickle_read)
    test_m = pickle.load(pickle_read_M)
    max_a_file = open("Max_Antigen.p",'r')
    min_a_file = open("Min_Antigen.p",'r')
    max_antigen = pickle.load(max_a_file)
    min_antigen = pickle.load(min_a_file)
    distancia_b = 1000
    distancia_m = 1000
    distancia = 0
    antigenos += antigenos2
    certo = 0
    errado = 0
    total = 0
    decisao = 'B' 
    for elemento in antigenos:
        for i in range(cord):
            elemento.coordenadas[i] = float(elemento.coordenadas[i])
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

    return str(float(certo)/float(total) * 100)+"%\n"

        

    ''' i = 1
    for element in test_pickle:
        print str(i) +" - "
        print element
        i += 1 '''
if __name__ == "__main__":
    START_TIME = time.time()
    filename = 'resultados.txt'
    try:
        f = open(filename,"r+")
    except IOError:
        f = open(filename,"w")
    f.seek(0,2)
    f.write('\n')

    for i in range(5):
        main_b(False)
        main_m(False)
        f.write(main())
    f.close()            
    print time.time()-START_TIME