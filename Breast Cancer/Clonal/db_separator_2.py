import time
import random
import copy
from math import floor

def ler_base(filename):
    f = open(filename, 'rU')
    total_data = []
    for line in f:
        data = line.split(',')
        data = data[1:12]
        total_data.append(data)
    return total_data
def data_shuffle(filename):
    """
    Retorna os dados randomizados lidos do dataset
    """
    new_data = ler_base(filename)
    random.shuffle(new_data)
    random.shuffle(new_data)
    random.shuffle(new_data)
    return new_data
def separate_data(data, filename, percentage = 0.5):
    '''
    Separa dados em Training e Test
    '''
    tam = len(data)
    tam = int(floor((tam*0.7)) + 1)
    data_training = data[:tam]
    data_test = data[tam:]
    return (data_training,data_test)
def separate_MB(data, filename, train = False):    
    if train:
        f_b = open(filename+'.outcome_B_training','w')
        f_m = open(filename+'.outcome_M_training','w')
        for element in data:
            if element[0] == 'B':
                phrase = ','.join(element)
                phrase += '\n'
                f_b.write(phrase)
            if element[0] == 'M':
                phrase = ','.join(element)
                phrase += '\n'
                f_m.write(phrase)
    else:
        f_b = open(filename+'.outcome_B','w')
        f_m = open(filename+'.outcome_M','w')
        for element in data:
            if element[0] == 'B':
                phrase = ','.join(element)
                phrase += '\n'
                f_b.write(phrase)
            if element[0] == 'M':
                phrase = ','.join(element)
                phrase += '\n'
                f_m.write(phrase)
    f_b.close()
    f_m.close()





if __name__ == '__main__':
    START_TIME = time.time()
    name = 'wdbc.data'
    n_data = data_shuffle(name)
    data_train,data_test = separate_data(n_data,name,0.7)
    separate_MB(data_train,name,train = True)
    separate_MB(data_test,name)
    print time.time()-START_TIME
