
from random import randint
import time
import math

class Antibody():
    def __init__(self, id=0, paratope=(0.5, 0.5), affinity=0):
        self.affinity = affinity
        self.relative_affinity = 1
        self.paratope = paratope
        self.id=id

    def set_relative_affinity(self, affinity):
        self.relative_affinity = affinity
    
    def __str__(self):
        s = 'Antibody: {}\n\tParatope: {}\n\tAffinity: {}\n'.format(self.id, self.paratope, self.affinity)
        return s

    def mutate(self, value, notify=False):
        '''Mutates paratope. Value is a percentage radius change'''
        x = self.paratope
        value = math.fabs(value)
        min, max = int(1000*(1-value)), int(1000*(1+value))     
        multiplyer = randint(min, max)/1000.0, randint(min, max)/1000.0 
        self.paratope = self.paratope[0]*multiplyer[0], self.paratope[1]*multiplyer[1] 
        self.update_affinity()
        if notify: print '\t{} * [{}~{}] = {} --> {} == {}'.format(x, min/1000.0, max/1000.0, multiplyer, self.paratope, self.affinity)
    
    def update_affinity(self):
        x, y = self.paratope
        pi = math.pi
        
        self.affinity = 0.26*(x*x + y*y) - 0.48*(x*y)

class Antigen():
    def __init__(self, id=0, epitope=None):
        self.epitope = epitope

def generate_antibodies(n, ag):
    antibodies = []
    for num in range(n):
        paratope = randint(-10000, 10000)/1000.0, randint(-10000, 10000)/1000.0 
        ab = Antibody(id=num, paratope=paratope)
        ab.update_affinity()
        antibodies.append(ab)
    return antibodies

def sort_key(antibody):
    return antibody.affinity

def main(notify=False):
    # constants
    ab_count = 200
    affinity_threshold = 0.2 # percentage
    ag = Antigen()

    antibodies = generate_antibodies(ab_count, ag)
    
    if notify:
        for antibody in antibodies:
            print antibody
        print '#'*30

    for iteration in range(10000):
        highest_aff_abs = []    
        antibodies = sorted(antibodies, key=sort_key, reverse=False)
        
        # select the highest affinity
        max_ab_affinity = antibodies[0].affinity
        for antibody in antibodies:
            if antibody.affinity >= max_ab_affinity * affinity_threshold:
                highest_aff_abs.append(antibody)
        
        # getting sum; will be used for weighs
        affs_sum = 0
        for antibody in highest_aff_abs:
            affs_sum += antibody.affinity

        # clones the best, up to ab_count
        cloned_highest_affs_abs = []
        for antibody in highest_aff_abs:
            n_clones = antibody.affinity*ab_count/affs_sum + 1
            for n in range(int(n_clones)):
                if len(cloned_highest_affs_abs) == ab_count:
                    break
                id = antibody.id
                aff = antibody.affinity
                paratope = antibody.paratope
                cloned_highest_affs_abs.append(Antibody(id=id, paratope=paratope, affinity=aff))

        # mutates clones and adds it to antibodies
        max_aff = highest_aff_abs[0].affinity
        for cloned_ab in cloned_highest_affs_abs:
            value = 1.01 - cloned_ab.affinity*1.0/max_aff
            cloned_ab.mutate(value, notify)
            antibodies.append(cloned_ab)

        antibodies = sorted(antibodies, key=sort_key, reverse=False)[:ab_count]      

    return antibodies[0]

if __name__ == '__main__':
    START_TIME = time.time()
    print main()
    print time.time()-START_TIME
