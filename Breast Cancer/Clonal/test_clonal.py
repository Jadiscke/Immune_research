import pickle
import csv
from ClonalgClassifier import Anticorpo, Antigeno, carregar, ler, normalizar
from ClonalgClassifier import distancia_euclidiana, main_b
from ClonalgClassifier_M import main_m
import time 
from db_separator_2 import db_separator
def voto(n_votos,at,test_b,test_m):
    #Retorna o chute do se antigenos se baseando nos votos
    #dos n melhores chutes
    vot_b = 0
    vot_m = 0
    distancias_tipo = []
    for anticorpo in test_b:
        distancia = distancia_euclidiana(at,anticorpo)
        distancias_tipo.append((distancia,'B'))    
    for anticorpo in test_m:
        distancia = distancia_euclidiana(at,anticorpo)
        distancias_tipo.append((distancia,'M')) 
    distancias_tipo.sort(key = lambda tup: tup[0])
    for i in range(n_votos):
        if distancias_tipo[i][1] == 'B':
            vot_b += 1
        else:
            vot_m += 1
    if vot_b > vot_m:
        decisao = 'B'
    else:
        decisao = 'M'
    return decisao

        
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
def tester(max_antigen,min_antigen,antigenos,test_n):
    cord = 10
    n = 3
    imp  = []
    #Carregar anticorpos
    filename = "Memory.p"
    filename2 = "Memory_M.p"
    pickle_read = open(filename,"r")
    pickle_read_M = open(filename2,"r")
    test_b = pickle.load(pickle_read) #anticorpos benignos
    test_m = pickle.load(pickle_read_M) #anticorpos malignos
   
    distancia = 0
    certo = 0
    errado = 0
    total = 0 
   
    

    
    for at in antigenos:
        decisao = maior_afinidade(at,test_b,test_m)
        print "Decisao - "+ decisao +" Original - " + at.tipo
        if decisao == at.tipo:
            certo += 1
        else:
            errado += 1
        total += 1
    #imprimir ="Maior afinidade - "+str(float(certo)/float(total) * 100)+"\t"
    percent = float(certo)/float(total) * 100
    #imp.append(imprimir)
    imp.append(percent)
    #Zerar
    certo = 0
    errado = 0
    total = 0 
    for at in antigenos:
        decisao = voto(n,at,test_b,test_m)
        if decisao == at.tipo:
            certo += 1
        else:
            errado += 1
        total += 1
    
    #imprimir2 = str(n)+" votos - "+str(float(certo)/float(total) * 100)+"\t"
    percent2 = float(certo)/float(total) * 100
    #imp.append(imprimir2)
    imp.append(percent2)
    #Zerar
    certo = 0
    errado = 0
    total = 0 
    n += 2
    for at in antigenos:
        decisao = voto(n,at,test_b,test_m)
        if decisao == at.tipo:
            certo += 1
        else:
            errado += 1
        total += 1
    

    #imprimir3 = str(n)+" votos - "+str(float(certo)/float(total) * 100)+"\t"
    percent3 = float(certo)/float(total) * 100
    #imp.append(imprimir3)
    imp.append(percent3)
    #Zerar
    certo = 0
    errado = 0
    total = 0 
    n += 2
    for at in antigenos:
        decisao = voto(n,at,test_b,test_m)
        if decisao == at.tipo:
            certo += 1
        else:
            errado += 1
        total += 1

    #imprimir4 = str(n)+" votos - "+str(float(certo)/float(total) * 100)+"\n"
    percent4 = float(certo)/float(total) * 100
    #imp.append(imprimir4)
    imp.append(percent4)
    return imp

        

    ''' i = 1
    for element in test_pickle:
        print str(i) +" - "
        print element
        i += 1 '''
if __name__ == "__main__":
    START_TIME = time.time()
    filename = 'resultados.csv'
    try:
        f = open(filename,"rb+")
    except IOError:
        f = open(filename,"wb")
    writer = csv.writer(f,delimiter = ' ',quotechar ='|', quoting = csv.QUOTE_MINIMAL )
    f.seek(0,2)
    i = 0
    n = 10
    results = [0.0,0.0,0.0,0.0]
    for i in range(n):
        db_separator()
        main_b(False)
        main_m(False)
        #Abrir Referencia de Antigenos
        max_a_file = open("Max_Antigen.p",'r')
        min_a_file = open("Min_Antigen.p",'r')
        #Carregar Referencia de antigenos para normalizacao
        max_antigen = pickle.load(max_a_file)
        min_antigen = pickle.load(min_a_file)
        max_a_file.close()
        min_a_file.close()
        #Carregar antigenos
        antigenos = carregar("wdbc.data.outcome_B")
        antigenos2 = carregar("wdbc.data.outcome_M")
        antigenos += antigenos2
        cord = 10
        # String -> Float
        k = 0
        for elemento in antigenos:
            for k in range(cord):
                elemento.coordenadas[k] = float(elemento.coordenadas[k])
        #Normalizando todos os antigenos
        k = 0
        for elemento in antigenos:
            for k in range(cord):
                elemento.coordenadas[k] = normalizar(elemento.coordenadas[k],min_antigen[k],max_antigen[k])

        antigeno = antigenos[:]
        c = 0
        r = tester(max_antigen,min_antigen,antigeno,i)
        writer.writerow(str(i)+"-"+str(r))
        for c in range(4):
            results[c] += r[c]
            c  += 1
        i+= 1
    i = 0
    writer.writerow("Media")
    for i in range(4):
        results[i] = results[i]/n
    writer.writerow(str(results))
    print results
    f.close()            
    print time.time()-START_TIME