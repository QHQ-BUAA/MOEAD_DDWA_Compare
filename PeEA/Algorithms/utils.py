import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
import openpyxl
import json
import configparser

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

def input_json(mk_code,cycle_code,EP):
    file_path = r'.\result'
    json_name = f'PeEA_data{mk_code}.json'
    json_path = os.path.join(file_path, json_name)
    if os.path.exists(json_path):
        with open(json_path, 'r') as json_file:
            existing_data = json.load(json_file)
    else:
        existing_data = {}
    arr2 = []
    for i in EP:
        temp_fit = i.function_value
        arr2.append(temp_fit)
    my_set = set(map(tuple, arr2))
    my_list = list(map(list, my_set))
    if cycle_code in existing_data:
        existing_data[cycle_code].extend(my_list)
    else:
        existing_data[cycle_code] = my_list
    with open(json_path, 'w') as json_file:
        json.dump(existing_data, json_file)

# 3-objective
def Tri_Dominate(Pop1,Pop2):
    '''
    :param Pop1:
    :param Pop2:
    :return: If Pop1 dominate Pop2, return True
    '''
    if (Pop1.function_value[0] < Pop2.function_value[0] and Pop1.function_value[1] < Pop2.function_value[1]) or \
        (Pop1.function_value[0] <= Pop2.function_value[0] and Pop1.function_value[1] < Pop2.function_value[1]) or \
        (Pop1.function_value[0] < Pop2.function_value[0] and Pop1.function_value[1] <= Pop2.function_value[1]) :
        return True
    else:
        return False

def object2_Dominate(list1, list2):
    '''
    :param list1:
    :param list2:
    :return: If Pop1 dominate Pop2, return True
    '''
    if (list1[0] < list2[0] and list1[1] < list2[1] or
        (list1[0] <= list2[0] and list1[1] < list2[1])or
        (list1[0] < list2[0] and list1[1] <= list2[1] )):
        return True
    else:
        return False
def dominate(num1,num2):
    if num1<num2:
        return True
    else:
        return False



# 3-objective
def TriPlot_NonDominatedSet(ax,color,EP,Instance,Algo_name,cpu_t):
    x = []
    y = []
    z=[]
    for i in range(len(EP)):
        x.append(EP[i].function_value[0])
        y.append(EP[i].function_value[1])
        z.append(EP[i].function_value[2])
    ax.scatter(x, y, z,color=color,label=str(Algo_name))
    ax.set_zlabel('Max Machine Load', fontdict={'size': 15, 'color': 'red'})
    ax.set_ylabel('Total Machine Load', fontdict={'size': 15, 'color': 'red'})
    ax.set_xlabel('Makespan', fontdict={'size': 15, 'color': 'red'})
    ax.view_init(elev=46, azim=33)
    plt.title('Instance: '+Instance+' '+'Algo_name: '+Algo_name+' '+'CPU(s): '+str(cpu_t))


def fast_non_dominated_sort(Pop):
    S=[[] for i in range(len(Pop))]
    front = [[]]
    n=[0 for i in range(len(Pop))]
    rank = [0 for i in range(len(Pop))]
    for p in range(len(Pop)):
        S[p]=[]
        n[p]=0
        for q in range(len(Pop)):
            if Tri_Dominate(Pop[p],Pop[q]):
                if q not in S[p]:
                    S[p].append(q)
            elif Tri_Dominate(Pop[q],Pop[p]):
                n[p] = n[p] + 1
        if n[p]==0:
            rank[p] = 0
            if p not in front[0]:
                front[0].append(p)
    i = 0
    while(front[i] != []):
        Q=[]
        for p in front[i]:
            for q in S[p]:
                n[q] =n[q] - 1
                if( n[q]==0):
                    rank[q]=i+1
                    if q not in Q:
                        Q.append(q)
        i = i+1
        front.append(Q)
    del front[len(front) - 1]
    NDSet=[]
    for Fi in front:
        NDSeti=[]
        for pi in Fi:
            NDSeti.append(Pop[pi])
        NDSet.append(NDSeti)
    return NDSet

def bi_VGM(Pop_size):
    delta=1/Pop_size
    w=[]
    w1=0
    while w1<=1:
        w2=1-w1
        w.append([w1,w2])
        w1+=delta
    return w


def get_inputInfo(key,Ihang,Jlie):
    config = configparser.ConfigParser()
    config.read(f'E:\PaperCode\Result\config\dauzere9-18\Ener\{Jlie}parameter{Ihang}.ini')

    ma_info = config.get(key, 'ma_dic')
    ma_dic = json.loads(ma_info)
    job_info = config.get(key, 'job_dic')
    job_dic = json.loads(job_info)
    return job_dic,ma_dic

def cal_fitness(PopObj):
    Zmin = np.min(PopObj, axis=0)
    M = PopObj.shape[1]
    W = np.zeros((M, M)) + 1e-6
    np.fill_diagonal(W, 1)
    N, M = PopObj.shape
    ASF = np.zeros((N, M))
    for i in range(M):
        ASF[:, i] = np.max((PopObj - Zmin) / W[i, :], axis=1)
    extreme = np.argmin(ASF, axis=0)  # 'extreme' is the extreme solutions
    Hyperplane = np.dot(np.linalg.pinv(PopObj[extreme, :]), np.ones(M))
    # Hyperplane = np.linalg.solve(PopObj[extreme, :], np.ones(M))  # Solves the linear equation
    a = 1 / Hyperplane
    PopObj = (PopObj - np.tile(Zmin, (N, 1))) / np.tile(a - Zmin, (N, 1))
    w_nad = a
    Zmin = np.min(PopObj, axis=0)
    ASF_q = np.zeros((N, M))
    for i in range(M):
        ASF_q[:, i] = np.max((PopObj - np.tile(Zmin, (N, 1))) / np.tile(w_nad, (N, 1)), axis=1)
    key_point = np.argmin(ASF_q, axis=0)
    key_point_indices = key_point.astype(int)
    selected_rows = PopObj[key_point_indices, :]
    q = np.sqrt(np.sum(np.square(selected_rows))) * np.sqrt(M)
    if q <= 1:
        fx = np.sum(PopObj - np.tile(Zmin, (N, 1)), axis=1)
    else:
        fx = np.max(PopObj - np.tile(Zmin, (N, 1)), axis=1)
    Distance = np.full((N, N), np.inf)
    for i in range(N):
        for j in range(N):
            if i != j:
                denominator = 1 + np.minimum(PopObj[i, :], PopObj[j, :])
                numerator = np.maximum(PopObj[i, :], PopObj[j, :]) - np.minimum(PopObj[i, :], PopObj[j, :])
                dis = numerator / denominator
                Distance[i, j] = np.sum(dis)
    Distance_sorted = np.sort(Distance, axis=1)
    k = int(np.floor(np.sqrt(N)))
    d = 1 / (Distance_sorted[:, k] + 2)
    Fitness = fx.flatten() + d
    return Fitness, extreme
