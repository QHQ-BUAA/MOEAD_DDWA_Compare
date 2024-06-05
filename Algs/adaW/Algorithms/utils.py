import os
import numpy as np
import pickle
import json

def Instance(file):
    '''
    :param file: Instance of FJSP
    example:
    n,m=3, 3 jobs, 3 machines
    PT=[[   # jobs 1
            [1,2],  # operation 1:  [1,2] indicatas for the first and second available machine's processing time
            [1,3]],
        [   # jobs 2
            [2,3],
            [3,2],]]
    MT=[[   # jobs 1
            [1,3], # operation 1:  [1,2] indicatas for the first and second available machine's index
            [3,2]],
        [   # jobs 2
            [1,2],
            [2,3]]]
    ni=[2,2] # the first job has 2 operations, the second job has 2 operations
    '''

    with open(file,"rb") as fb:
        I=pickle.load(fb)
    if type(I) == str:
        I = json.loads(I)
    n,m,PT,MT,ni=I['n'],I['m'],I['processing_time'],I['Processing machine'],I['Jobs_Onum']

    return n,m,PT,MT,ni

def Tchebycheff(x, z, lambd):
    '''
    :param x: Popi
    :param z: the reference point
    :param lambd: a weight vector
    :return: Tchebycheff objective
    '''
    Gte = []
    for i in range(len(x.fitness)):
        Gte.append(np.abs(x.fitness[i]-z[i]) * lambd[i])#
    return np.max(Gte)


def input_json(mk_code,cycle_code,EP):
    file_path = r'.\result'
    json_name = f'IGD_data{mk_code}.json'
    json_path = os.path.join(file_path, json_name)
    if os.path.exists(json_path):
        with open(json_path, 'r') as json_file:
            existing_data = json.load(json_file)
    else:
        existing_data = {}
    arr2 = []
    for i in EP:
        temp_fit = i.fitness
        arr2.append(temp_fit)
    my_set = set(map(tuple, arr2))
    my_list = list(map(list, my_set))
    if cycle_code in existing_data:
        existing_data[cycle_code].extend(my_list)
    else:
        existing_data[cycle_code] = my_list
    with open(json_path, 'w') as json_file:
        json.dump(existing_data, json_file)

def Neighbor(lambd, T):
    B = []
    for i in range(len(lambd)):
        temp = []
        for j in range(len(lambd)):
            distance = np.sqrt((lambd[i][0] - lambd[j][0])**2 + (lambd[i][1] - lambd[j][1])**2 )
            temp.append(distance)
        res = np.argsort(temp)
        B.append(res[:T])
    return B


def Tri_Dominate(Pop1,Pop2):
    '''
    :param Pop1:
    :param Pop2:
    :return: If Pop1 dominate Pop2, return True
    '''
    if (Pop1.fitness[0]<Pop2.fitness[0] and Pop1.fitness[1]<Pop2.fitness[1]) or \
        (Pop1.fitness[0] <= Pop2.fitness[0] and Pop1.fitness[1] < Pop2.fitness[1]) or \
        (Pop1.fitness[0] < Pop2.fitness[0] and Pop1.fitness[1] <= Pop2.fitness[1]) :
        return True
    else:
        return False

def bi_VGM(Pop_size):
    delta=1/Pop_size
    w=[]
    w1=0
    while w1<=1:
        w2=1-w1
        w.append([w1,w2])
        w1+=delta
    return w



