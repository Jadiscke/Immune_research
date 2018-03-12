import pickle
from ClonalgClassifier import Anticorpo, Antigeno, carregar, ler, normalizar
from ClonalgClassifier import distancia_euclidiana, main_b
from ClonalgClassifier_M import main_m
import time 
from db_separator_2 import db_separator

def voto(n_votos,at,test_b,test_m):
    '''
    Retorna o chute do se antigenos se baseando nos votos
    dos n melhores chutes
    '''
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
    '''
    retorna o chute da antigeno se baseando na celula de memoria mais proxima
    '''
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
def tester(max_antigen,min_antigen,antigenos):
    n = 3

    #Carregar anticorpos
    filename = "Memory.p"
    filename2 = "Memory_M.p"
    pickle_read = open(filename,"r")
    pickle_read_M = open(filename2,"r")
    test_b = pickle.load(pickle_read) #anticorpos benignos
    test_m = pickle.load(pickle_read_M) #anticorpos malignos
   
    vp = 0.0 #Verdadeiro Positivo
    vn = 0.0 #Verdadeiro Negativo
    fp = 0.0 #Falso Positivo
    fn = 0.0 #Falso Negativo   
    certo = 0.0
    errado = 0.0
    total = 0.0

    
    #Teste for maior afinidade
    for at in antigenos:
        decisao = maior_afinidade(at,test_b,test_m)
        if decisao == at.tipo:
            if decisao == 'B':
                vp += 1
            else:
                vn += 1
            certo += 1
            
        else:
            if decisao == 'B':
                fp += 1
            else:
                fn += 1
            errado += 1
        total += 1
    precision = (vp)/(vp + fp) 
    recall = (vp)/(vp + fn)
    f1score = 2 *((precision * recall)/(precision + recall))
    percentual = certo/total
    imprimir = [percentual,precision,recall,f1score]
    #imprimir ="Per- "+str("%.2f" % (percentual*100))+"%\tPrec- "+str("%.2f" %(precision*100)) + "%\tRec- "+str("%.2f" % (recall*100))+ "%\tF1- "+str("%.2f" % (f1score*100))+"%\n"
    #Zerar
    vp = 0.0 
    vn = 0.0
    fp = 0.0
    fn = 0.0
    certo = 0.0
    errado = 0.0
    total = 0.0 
    for at in antigenos:
        decisao = voto(n,at,test_b,test_m)
        if decisao == at.tipo:
            if decisao == 'B':
                vp += 1
            else:
                vn += 1
            certo += 1
        else:
            if decisao == 'B':
                fp += 1
            else:
                fn += 1
            errado += 1
        total += 1
    
    precision = (vp)/(vp + fp)
    recall = (vp)/(vp + fn)
    f1score = 2 *((precision * recall)/(precision + recall))
    percentual = certo/total
    imprimir2 = [percentual,precision,recall,f1score]
    #imprimir2 ="Per- "+str("%.2f" % (percentual*100))+"%\tPrec- "+str("%.2f" %(precision*100)) + "%\tRec- "+str("%.2f" % (recall*100))+ "%\tF1- "+str("%.2f" % (f1score*100))+"%\n"
    #Zerar
    vp = 0.0 
    vn = 0.0
    fp = 0.0
    fn = 0.0
    certo = 0.0
    errado = 0.0
    total = 0.0 
    n += 2
    for at in antigenos:
        decisao = voto(n,at,test_b,test_m)
        if decisao == at.tipo:
            if decisao == 'B':
                vp += 1
            else:
                vn += 1

            certo += 1
        else:
            if decisao == 'B':
                fp += 1
            else:
                fn += 1

            errado += 1
        total += 1
    
    precision = (vp)/(vp + fp)
    recall = (vp)/(vp + fn)
    f1score = 2 *((precision * recall)/(precision + recall))
    percentual = certo/total
    imprimir3 = [percentual,precision,recall,f1score]
    #imprimir3 ="Per- "+str("%.2f" % (percentual*100))+"%\tPrec- "+str("%.2f" %(precision*100)) + "%\tRec- "+str("%.2f" % (recall*100))+ "%\tF1- "+str("%.2f" % (f1score*100))+"%\n"
    
    #Zerar
    vp = 0.0 
    vn = 0.0
    fp = 0.0
    fn = 0.0
    certo = 0.0
    errado = 0.0
    total = 0.0 
    n += 2
    for at in antigenos:
        decisao = voto(n,at,test_b,test_m)
        if decisao == at.tipo:
            if decisao == 'B':
                vp += 1
            else:
                vn += 1
            certo += 1
        else:
            if decisao == 'B':
                fp += 1
            else:
                fn += 1

            errado += 1
        total += 1
    
    precision = (vp)/(vp + fp)
    recall = (vp)/(vp + fn)
    f1score = 2 *((precision * recall)/(precision + recall))
    percentual = certo/total
    imprimir4 = [percentual,precision,recall,f1score]
    #imprimir4 ="Per- "+str("%.2f" % (percentual*100))+"%\tPrec- "+str("%.2f" %(precision*100)) + "%\tRec- "+str("%.2f" % (recall*100))+ "%\tF1- "+str("%.2f" % (f1score*100))+"%\n"
    
    return (imprimir, imprimir2,imprimir3,imprimir4)

        

if __name__ == "__main__":
    START_TIME = time.time()
    #Arquivos para salvar dados
    #filename1,filename2,filename3,filename4 = 'Afinidade.txt','3votos.txt','5votos.txt','7votos.txt'
    
    # try:
    #     f1,f2,f3,f4 = open(filename1,"r+"),open(filename2,"r+"),open(filename3,"r+"),open(filename4,"r+")
    # except IOError:
    #     f1,f2,f3,f4 = open(filename1,"w"),open(filename2,"w"),open(filename3,"w"),open(filename4,"w")
    # #Escrever no final do Arquivo
    # f1.seek(0,2)
    # f2.seek(0,2)
    # f3.seek(0,2)
    # f4.seek(0,2)
    for i in range(10):
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
        for elemento in antigenos:
            for k in range(cord):
                elemento.coordenadas[k] = float(elemento.coordenadas[k])
        #Normalizando todos os antigenos
        for elemento in antigenos:
            for k in range(cord):
                elemento.coordenadas[k] = normalizar(elemento.coordenadas[k],min_antigen[k],max_antigen[k])

    
        #Testando e escrevendo no arquivo
        t1,t2,t3,t4 = tester(max_antigen,min_antigen,antigenos)
        print i
        print t1
        print t2
        print t3
        print t4
        print 20 * '-'
        #f1.write(t1)
        #f2.write(t2)
        #f3.write(t3)
        #f4.write(t4)
    #f1.close()
    #f2.close()  
    #f3.close()  
    #f4.close()              
    print time.time()-START_TIME