#utf-8
import pdb
import random
import copy
from Algorithms.utils import *
from Algorithms.Popi import *
import numpy as np
from scipy.spatial.distance import pdist, squareform

from PeEA.Algorithms.utils import fast_non_dominated_sort


class Algorithms:
    def __init__(self,args,job_dic,ma_dic):
        self.means_m =args.means_m
        self.job_dic = job_dic
        self.ma_dic = ma_dic
        self.args =args
        self.N_elite =args.N_elite  #精英数量
        self.Pop_size =args.pop_size
        self.gene_size =args.gene_size
        self.pc_max =args.pc_max
        self.pm_max =args.pm_max
        self.pc_min = args.pc_min
        self.pm_min = args.pm_min
        self.pc = args.pc_max
        self.pm = args.pm_max
        self. T =args.T
        if self.means_m >1:
            self.p_GS =args.p_GS
            self.p_LS = args.p_LS
            self.p_RS = args.p_RS

        if self.means_m >1:
            self.Chromo_setup()  #生成了Ms_List,Os_list,half_len_chromo(工序数)，J_site(索引)

        self.Best_JS =None
        self.Best_Cmax =1e+20
        self.C_end =[]
        self.Pop = []       # Population
        # self._lambda =bi_VGM(self.Pop_size)     # bi-objective weight vectors
        self._lambda = bi_VGM(self.args.H)    #生成权重向量  # Tri-objective weight vectors   # H：沿每个目标坐标考虑的划分数量
        self._z =[]          # the reference point
        self.min_mt = 0
        self.mint_pop = []


    def Chromo_setup(self):
        self.os_list = []
        for i in range(len(self.args.O_num)):
            self.os_list.extend([i for _ in range(self.args.O_num[i])])
        self.half_len_chromo =len(self.os_list)
        self.ms_list =[] #每个工序有几个机器可以加工
        self.J_site =[]  # 方便后面的位置查找
        for i in range(len(self.args.Processing_Machine)):
            for j in range(len(self.args.Processing_Machine[i])):
                self.ms_list.append(len(self.args.Processing_Machine[i][j]))
                self.J_site.append((i,j))

    def random_initial(self):
        for i in range(self.Pop_size):  #随机初始化概率 * 种群数  #在这里可以生成想要多少随机的染色体
            Pop_i =[]
            random.shuffle(self.os_list)
            Pop_i.extend(copy.copy(self.os_list))
            ms =[]
            for i in self.ms_list:
                ms.append(random.randint(0, i -1))#生成一个在范围 [a,b] 内的随机整数。其中，a是范围的下限，b是范围的上限。这个函数的返回值是一个整数，它的值在 [a,b] 范围内，包括a和b。
            Pop_i.extend(ms)  #相当于生成了一个随机初始化的染色体
            Pop_i =Popi(self.args,Pop_i,self.job_dic, self.ma_dic,self.J_site,self.half_len_chromo)  #适应度返回的分别是 最大制造跨度   总负荷  某台机器的最大负荷
            if self._z==[]:
                self._z =Pop_i.function_value
            else:
                for j in range(2):
                    if self._z[j ] >Pop_i.function_value[j]:
                        self._z[j ] =Pop_i.function_value[j]  #这里有个bug，每次都会更新Pop[0]的fitness，无法更改
            self.Pop.append(Pop_i)

    def MatingSelection(self):
        n = len(self.Pop)
        funva_list = []
        for ni in range(n):
            funva_list.append(self.Pop[ni].function_value)
        fitness,extreme = cal_fitness(np.array(funva_list))
        for ni in range(n):
            self.Pop[ni].fitness = fitness[ni]
        p1 = random.choices(range(n), k=n)  #随机选两个父代
        p2 = random.choices(range(n), k=n)
        for pi in range(n): #else:如果1、2互不支配，则比较fitness
            if Tri_Dominate(self.Pop[p1[pi]],self.Pop[p2[pi]]):
                self.Pop[p2[pi]] = self.Pop[p1[pi]]
            elif Tri_Dominate(self.Pop[p2[pi]],self.Pop[p1[pi]]):
                self.Pop[p1[pi]] = self.Pop[p2[pi]]
            else:
                if dominate(self.Pop[p1[pi]].fitness,self.Pop[p2[pi]].fitness):
                    self.Pop[p2[pi]] = self.Pop[p1[pi]]
                elif dominate(self.Pop[p2[pi]].fitness,self.Pop[p1[pi]].fitness):
                    self.Pop[p1[pi]] = self.Pop[p2[pi]]

    def optimization(self,R_pop):
        n = len(R_pop)
        funva_list = []
        for ni in range(n):
            funva_list.append(R_pop[ni].function_value)
        funva_array = np.array(funva_list)
        fitness, extreme = cal_fitness(funva_array)
        for ni in range(n):
            R_pop[ni].fitness = fitness[ni]
        extreme = np.unique(extreme)
        Angle = np.arccos(1 - squareform(pdist(funva_array, 'cosine')))
        np.fill_diagonal(Angle, np.inf)
        Remain = list(range(len(funva_array)))
        for e in extreme:
            if e in Remain:
                Remain.remove(e)
        while len(Remain) > (self.Pop_size - len(extreme)):
            AngleSubset = Angle[Remain, :][:, Remain]
            sortA = np.sort(AngleSubset, axis=1)
            rank1 = np.argsort(AngleSubset, axis=1)
            rank2 = np.argsort(sortA, axis=0)[0, :]
            A = rank2[0]
            B = rank1[A, 0]

            if fitness[Remain[A]] > fitness[Remain[B]]:
                del Remain[A]
            else:
                del Remain[B]
        indices = Remain + list(extreme)
        R_pop = [R_pop[i] for i in indices]
        return R_pop
    '''
    工序交叉：
    （1）POX
    （2）Job_based_Crossover
    '''

    # POX:precedence preserving order-based crossover
    def POX(self,p1, p2):
        jobsRange = range(0, self.args.n)
        sizeJobset1 = random.randint(1, self.args.n)
        jobset1 = random.sample(jobsRange, sizeJobset1)
        o1 = []
        p1kept = []
        for i in range(len(p1)):
            e = p1[i]
            if e in jobset1:
                o1.append(e)
            else:
                o1.append(-1)
                p1kept.append(e)
        o2 = []
        p2kept = []
        for i in range(len(p2)):
            e = p2[i]
            if e in jobset1:
                o2.append(e)
            else:
                o2.append(-1)
                p2kept.append(e)
        for i in range(len(o1)):
            if o1[i] == -1:
                o1[i] = p2kept.pop(0)
        for i in range(len(o2)):
            if o2[i] == -1:
                o2[i] = p1kept.pop(0)
        return o1, o2

    def Job_Crossover(self,p1,p2):
        jobsRange = range(0, self.args.n)
        sizeJobset1 = random.randint(0, self.args.n)
        jobset1 = random.sample(jobsRange, sizeJobset1)
        jobset2 = [item for item in jobsRange if item not in jobset1]
        o1 = []
        p1kept = []
        for i in range(len(p1)):
            e = p1[i]
            if e in jobset1:
                o1.append(e)
                p1kept.append(e)
            else:
                o1.append(-1)
        o2 = []
        p2kept = []
        for i in range(len(p2)):
            e = p2[i]
            if e in jobset2:
                o2.append(e)
                p2kept.append(e)
            else:
                o2.append(-1)
        for i in range(len(o1)):
            if o1[i] == -1:
                o1[i] = p2kept.pop(0)
        for i in range(len(o2)):
            if o2[i] == -1:
                o2[i] = p1kept.pop(0)
        return o1,o2

    '''
    工序变异：
    swap_mutation
    NB_mutation: neigborhood mutation
    '''
    def swap_mutation(self,p):
        pos1 = random.randint(0, len(p) - 1)
        pos2 = random.randint(0, len(p) - 1)
        if pos1 == pos2:
            return p
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        offspring = p[:pos1] + [p[pos2]] + \
                    p[pos1 + 1:pos2] + [p[pos1]] + \
                    p[pos2 + 1:]
        return offspring

    def MB_mutation(self,p1):
        D = len(p1)
        c1 = p1.copy()
        r = np.random.uniform(size=D)
        for idx1, val in enumerate(p1):
            if r[idx1] <= self.pm:
                idx2 = np.random.choice(np.delete(np.arange(D), idx1))
                c1[idx1], c1[idx2] = c1[idx2], c1[idx1]
        return c1

    def Crossover_Machine(self, CHS1, CHS2):
        T_r = [j for j in range(self.half_len_chromo)]
        r = random.randint(1, self.half_len_chromo)  # 在区间[1,T0]内产生一个整数r
        random.shuffle(T_r)
        R = T_r[0:r]  # 按照随机数r产生r个互不相等的整数
        # 将父代的染色体复制到子代中去，保持他们的顺序和位置
        for i in R:
            K, K_2 = CHS1[i], CHS2[i]
            CHS1[i], CHS2[i] = K_2, K
        return CHS1, CHS2

    def Mutation_Machine(self,CHS):
        T_r = [j for j in range(self.half_len_chromo)]
        r = random.randint(1, self.half_len_chromo)  # 在区间[1,T0]内产生一个整数r
        random.shuffle(T_r)
        R = T_r[0:r]  # 按照随机数r产生r个互不相等的整数
        for i in R:
            O_site =self.J_site[i]
            pt =self.args.Processing_Time[O_site[0]][O_site[1]]
            pt_find =pt[0]
            len_pt =len(pt ) -1
            k, m =1,0
            while k< len_pt:
                if pt_find > pt[k]:
                    pt_find = pt[k]
                    m = k
                k += 1
            CHS[i] = m
        return CHS

    def operator_NoFlexible(self, chs1, chs2):
        p1, p2 = chs1.CHS, chs2.CHS
        if random.random() < self.pc:  # wether crossover
            if random.random() < 0.5:
                p1, p2 = self.POX(p1, p2)
            else:
                p1, p2 = self.Job_Crossover(p1, p2)
        if random.random() < self.pm:  # wether chs1 mutation
            if random.random() < 0.5:
                p1 = self.swap_mutation(p1)
            else:
                p1 = self.MB_mutation(p1)
        if random.random() < self.pm:  # wether chs2 mutation
            if random.random() < 0.5:
                p2 = self.swap_mutation(p2)
            else:
                p2 = self.MB_mutation(p2)
        P1, P2 = Popi(self.args, p1,self.job_dic, self.ma_dic, self.J_site, self.half_len_chromo), \
                    Popi(self.args, p2,self.job_dic, self.ma_dic, self.J_site,self.half_len_chromo)
        return P1, P2

    def operator_Flexible(self, chs1, chs2):
        p1, p2 = chs1.CHS, chs2.CHS
        if random.random() < self.pc:  # wether crossover
            if random.random() < 0.5:
                p11, p21 = self.POX(p1[0:self.half_len_chromo], p2[0:self.half_len_chromo])
            else:
                p11, p21 = self.Job_Crossover(p1[0:self.half_len_chromo], p2[0:self.half_len_chromo])
            p12, p22 = self.Crossover_Machine(p1[self.half_len_chromo:], p2[self.half_len_chromo:])
            p11.extend(p12)
            p1 = p11
            p21.extend(p22)
            p2 = p21
        if random.random() < self.pm:  # wether chs1 mutation
            if random.random() < 0.5:
                p11 = self.swap_mutation(p1[0:self.half_len_chromo])
            else:
                p11 = self.MB_mutation(p1[0:self.half_len_chromo])
            p12 = self.Mutation_Machine(p1[self.half_len_chromo:])
            p11.extend(p12)
            p1 = p11
        if random.random() < self.pm:  # wether chs2 mutation
            if random.random() < 0.5:
                p21 = self.swap_mutation(p2[0:self.half_len_chromo])
            else:
                p21 = self.MB_mutation(p2[0:self.half_len_chromo])
            p22 = self.Mutation_Machine(p2[self.half_len_chromo:])
            p21.extend(p22)
            p2 = p21
        P1, P2 = Popi(self.args, p1,self.job_dic, self.ma_dic, self.J_site, self.half_len_chromo), \
            Popi(self.args, p2,self.job_dic, self.ma_dic, self.J_site, self.half_len_chromo)
        return P1, P2

    def offspring_Population(self):
        new_pop=[]
        while len(new_pop)<self.Pop_size:
            pop1,pop2=random.choice(self.Pop),random.choice(self.Pop)
            if self.means_m>1:
                new_pop1,new_pop2=self.operator_Flexible(pop1,pop2)
            else:
                new_pop1, new_pop2 = self.operator_NoFlexible(pop1, pop2)
            new_pop.extend([new_pop1,new_pop2])
        return new_pop


    '''
    NSGA-Ⅱ原理，讲的十分清楚：
    https://www.bilibili.com/video/BV1v34y1y728/?spm_id_from=333.999.top_right_bar_window_history.content.click&vd_source=3b20d8b54b0e354a9780ff9a740b3a1c
    '''
    def PeEA_main(self, Ihang, Jlie):
        # ----------------------------------Initialization---------------------------------
        self.Pop=[]
        if self.means_m > 1:  # Used to determine whether the machine is flexible
            self.random_initial()  # random Initial
        # ----------------------------------Iteration---------------------------------
        for i in range(self.gene_size):
            self.MatingSelection() #形成新的Pop，返回极端解
            new_pop=self.offspring_Population()
            R_pop=self.Pop+new_pop
            self.Pop =  self.optimization(R_pop)
        EP=fast_non_dominated_sort(self.Pop)[0]
        return EP
