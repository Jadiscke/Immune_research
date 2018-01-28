import time
import random
import copy

def read_data(filename):
    f = open(filename, 'rU')
    total_data = []
    for line in f:
        data = line.split(',')
        data[-1] = data[-1].strip()
        total_data.append(data)
    return total_data

def separate_on_outcome(filename, outcome_index, separate_training_data=False, training_percentage=0.5):
    data = read_data(filename)
    outcomes = {}
    for line in data:
        if line[outcome_index] in outcomes:
            outcomes[line[outcome_index]].append(line)
        else:
            outcomes[line[outcome_index]] = [line]
    
    if separate_training_data:
        outcomes_copy = copy.deepcopy(outcomes)
        for outcome in outcomes_copy:
            outcomes[outcome+'_training'] = []
            while len(outcomes[outcome+'_training'])*1.0/len(outcomes_copy[outcome]) < training_percentage:
                line = outcomes[outcome][random.randint(0, len(outcomes[outcome]) - 1)]
                if line not in outcomes[outcome+'_training']:
                    outcomes[outcome+'_training'].append(line)
                    outcomes[outcome].remove(line)
                

    return outcomes

def save_outcomes(outcomes, filename):
    for outcome in outcomes:
        f = open(filename+'.outcome_'+str(outcome), 'w')
        for line in outcomes[outcome]:
            line_str = ''
            for parameter in line:
                line_str += str(parameter) + ','
            line_str = line_str[:-1] + '\n'
            f.write(line_str)


if __name__ == '__main__':
    START_TIME = time.time()
    filename = 'wdbc.data'
    outcomes = separate_on_outcome(filename, 1, separate_training_data=True)
    save_outcomes(outcomes, filename)
    print time.time()-START_TIME
