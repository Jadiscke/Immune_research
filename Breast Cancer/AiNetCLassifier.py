from random import randint, uniform
import time
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import random
#Anticorpo
class Antibody():
    def __init__(self, id=0, paratope=[], affinity=0, neighbors=4):
        self.affinity = affinity
        self.relative_affinity = 1
        self.paratope = paratope
        self.id=id
        self.neighbors = neighbors

        for index, parameter in enumerate(paratope):
            if parameter == 'M': parameter = 1
            if parameter == 'B': parameter = 0
            setattr(self, 'parameter_'+str(index), float(parameter))
    #Print
    def __str__(self):
        s = 'Antibody: {}\n\tParatope: {}\n'.format(self.id, self.paratope)
        for index, attr in enumerate(dir(self)):
            if 'parameter_' not in attr:
                continue
            s += '\t{}: {}\n'.format(attr, getattr(self, attr))
        s += '\tAffinity: {}\n'.format(self.affinity)
        return s

    def k_neighbors_affinity(self, k, antigens, notify=False):
        START_TIME = time.time()

        ags_distances = {}
        for antigen in antigens:
            ags_distances[antigen.parameter_0] = euclidian_distance(self, antigen)

        ags_distances = sorted(ags_distances.items(), key=lambda x: x[1])[:k]
        neares_ags_ids = [item[0] for item in ags_distances]
        near_antigens = [antigen for antigen in antigens if antigen.parameter_0 in neares_ags_ids]
        
        parameters_sum = {}
        for antigen in near_antigens:
            for attr in dir(antigen):
                if 'parameter_' in attr and attr != 'parameter_0' and attr != 'parameter_1':
                    if attr not in parameters_sum:
                        parameters_sum[attr] = getattr(antigen, attr)
                    else:
                        parameters_sum[attr] += getattr(antigen, attr)
        
        for parameter in parameters_sum:
            parameters_sum[parameter] = parameters_sum[parameter]/k

        if notify:
            print 'K_RUNTIME:', time.time() - START_TIME

        self.affinity = 100.0/(1 + euclidian_distance(self, parameters_sum))

    def mutate(self, value, notify=False):
        for attr in dir(self):
            if 'parameter_' in attr and attr != 'parameter_0' and attr != 'parameter_1':
                rand = np.random.rand()* random.randint(-1, 1)
                while rand == 0:
                    rand = np.random.rand()* random.randint(-1, 1)
                setattr(self, attr, getattr(self, attr) * (1 + (value * rand)))

    def update_affinity(self):
        pass
#Antigeno    
class Antigen():
    def __init__(self, *args):
        for index, arg in enumerate(*args):
            if arg == 'M': arg = 1
            if arg == 'B': arg = 0
            setattr(self, 'parameter_' + str(index), float(arg))

    def __str__(self):
        s = 'Antigen {} Type {}:\n'.format(self.parameter_0, self.parameter_1)
        for index, attr in enumerate(dir(self)):
            if '__' in attr or attr == 'parameter_0' or attr == 'parameter_1':
                continue
            s += '\t{}: {}\n'.format(attr, getattr(self, attr))
        return s


def euclidian_distance(a, b):
    distance_sum = 0
    if type(b) is dict:
        for attr in dir(a):
            if 'parameter_' in attr and attr != 'parameter_0' and attr != 'parameter_1':
                distance_sum += (getattr(a, attr) - b[attr]) ** 2
    else:
        for attr in dir(a):
            if 'parameter_' in attr and attr != 'parameter_0' and attr != 'parameter_1':
                distance_sum += (getattr(a, attr) - getattr(b, attr)) ** 2
    return distance_sum ** 0.5

def read_data(filename):
    f = open(filename, 'rU')
    total_data = []
    for line in f:
        data = line.split(',')
        data[-1] = data[-1].strip()
        total_data.append(data)
    return total_data

def load_data(filename):
    data = read_data(filename)
    antigens = []
    for line in data:
        antigens.append(Antigen(line))
    return antigens

def generate_antibodies(n, ags, neighbors, update_affinity=True):
    antibodies = []
    for index, ag in enumerate(ags):
        paratope = []
        for attr in sorted(dir(ag)):
            if 'parameter_' in attr:
                paratope.append(getattr(ag, attr))
        ab = Antibody(index, paratope, 0, neighbors)
        if update_affinity:
            ab.k_neighbors_affinity(neighbors, ags)
        antibodies.append(ab)        
    return antibodies

def sort_key(antibody):
    return antibody.affinity

# o que seriam multiplos antigenos?
# qual funcao eu utilizo pra definir o numero de clones?
def main(notify=False):
    # constants
    decay_rate = -10
    neighbors = 5
    max_clones = 40
    ab_count = 20
    clones_count = 4
    affinity_threshold = 0.5 # percentage that gets OKed
    
    antigens = load_data('wdbc.data')
    init_antibodies = generate_antibodies(ab_count, antigens, neighbors, False)
    antibodies = []

    START_TIME = time.time()

    # initiate antibodies as random cases of ags
    if len(init_antibodies) > ab_count:
        for i in range(ab_count):
            ab = init_antibodies[randint(0, len(init_antibodies) - 1)]
            while ab in antibodies:
                ab = init_antibodies[randint(0, len(init_antibodies) - 1)]
            antibodies.append(ab)
        
    # gets initial affinities
    for antibody in antibodies:
        antibody.k_neighbors_affinity(neighbors, antigens)

    if notify:
        print '#'*20
        print 'INITIAL ANTIBODIES'
        for antibody in antibodies:
            print antibody
        print '#'*20

    for iteration in range(30): 
        print '\nITERATION: ', iteration
        antibodies = sorted(antibodies, key=sort_key, reverse=True)
        highest_aff_abs = []
        # Select the highest affinity
        max_ab_affinity = antibodies[0].affinity
        for antibody in antibodies:
            highest_aff_abs.append(antibody)
        
        # Getting sum; will be used for weighs
        affs_sum = 0
        for antibody in highest_aff_abs:
            affs_sum += antibody.affinity

        # Clones the best
        max_ab_affinity = highest_aff_abs[0].affinity
        min_ab_affinity = highest_aff_abs[-1].affinity
        for antibody in highest_aff_abs:
            normalized_affinity = 1
            if (max_ab_affinity-min_ab_affinity) > 0:
                normalized_affinity = (antibody.affinity-min_ab_affinity)*1.0/(max_ab_affinity-min_ab_affinity)
            n_clones = normalized_affinity*clones_count + 1
            for n in range(int(n_clones)):      
                # Creates child and mutates it.
                # If its better than the parent, go-ahead. Else, dump it
                cloned_ab = Antibody(id=antibody.id, paratope=antibody.paratope, affinity=antibody.affinity)
                value = (math.expm1(normalized_affinity*decay_rate) + 1)
                cloned_ab.mutate(value)
                cloned_ab.k_neighbors_affinity(neighbors, antigens)
                if cloned_ab.affinity >= antibody.affinity:
                    print "ITERATION", iteration, "THERE'S A BETTER CLONE: {} VERSUS {} MUTATION VALUE {}".format(antibody.affinity, cloned_ab.affinity, value)
                    antibodies.append(cloned_ab)
                else:
                    print "ITERATION", iteration, "WORSE CLONE: {} VERSUS {} MUTATION VALUE {}".format(antibody.affinity, cloned_ab.affinity, value)


        antibodies = sorted(antibodies, key=sort_key, reverse=True)
        antibodies = antibodies[:max_clones]

        for antibody in antibodies:
            print antibody.affinity,

        # removes antibodies that are too close
        antibodies_copy = antibodies[:]
        for antibody in antibodies_copy:
            for antibody2 in antibodies_copy:
                if euclidian_distance(antibody, antibody2) < 10:
                    if antibody.affinity < antibody2.affinity and antibody in antibodies:
                        antibodies.remove(antibody)


        if iteration < 100:
            # Adds new batch of randomized abs
            antibodies = antibodies + generate_antibodies(ab_count, antigens, neighbors)

        antibodies = sorted(antibodies, key=sort_key, reverse=True) 

    if notify:
        print '#'*20
        print 'FINAL ANTIBODIES'
        for antibody in antibodies:
            print antibody
        print '#'*20

    return antibodies[0]

if __name__ == '__main__':
    START_TIME = time.time()
    main()
    print time.time()-START_TIME
