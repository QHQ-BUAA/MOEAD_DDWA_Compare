#utf-8
import copy
import random
from adaW.Algorithms.utils import *
from adaW.Algorithms.Popi import *
from scipy.spatial.distance import cdist


class Algorithms:
    def __init__(self ,args,Ihang,JLie):
        self.Ihang = Ihang
        self.JLie = JLie
        self.means_m =args.means_m
        self.args =args
        self.N_elite =args.N_elite
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
            self.p_GS = args.p_GS
            self.p_LS = args.p_LS
            self.p_RS = args.p_RS
        else:
            self.p_GS = 0
            self.p_LS = 0
            self.p_RS = 1
        if self.means_m >1:
            self.Chromo_setup()
        else:
            self.Chromo_setup_0()
        self.Best_JS =None
        self.Best_Cmax =1e+20
        self.C_end =[]
        self.Pop = []       # Population
        self._lambda = bi_VGM(self.args.H)
        self._z =[]          # the reference point
        self.min_1fit = 0
        self.min_pop = []

    def weight_update(self,EP):
        populationFitness = []
        for pfi in range(len(self.Pop)):
            populationFitness.append(self.Pop[pfi].fitness)
        population_fitness = np.array(populationFitness)
        EPFitness = []
        for epi in range(len(EP)):
            EPFitness.append(EP[epi].fitness)
        my_set = set(map(tuple, EPFitness))
        my_list = list(map(list, my_set))
        EP_fitness = np.array(my_list)
        W = np.array(self._lambda)
        N = population_fitness.shape[0]
        all_fitness = np.vstack((population_fitness, EP_fitness))
        fmin = np.min(all_fitness, axis=0)
        fmax = np.max(all_fitness, axis=0)

        normalized_EP_fitness = (EP_fitness - fmin) / (fmax - fmin)
        normalized_population_fitness = (population_fitness - fmin) / (fmax - fmin)
        dis1 = cdist(normalized_EP_fitness, normalized_population_fitness)
        dis1_sorted = np.sort(dis1, axis=1)
        if normalized_EP_fitness.shape[0] > 1:
            dis2 = cdist(normalized_EP_fitness, normalized_EP_fitness)
            np.fill_diagonal(dis2, np.inf)
            dis2_sorted = np.sort(dis2, axis=1)
            niche_size = np.median(dis2_sorted[:, 1])
        else:
            niche_size = 1
        undeveloped_indices = np.where(dis1_sorted[:, 0] >= niche_size)[0]
        z = np.array(self._z)
        for undeveloped_idx in undeveloped_indices:
            undeveloped_fitness = EP_fitness[undeveloped_idx]
            new_weight = (undeveloped_fitness - z) / (np.sum(undeveloped_fitness) - np.sum(z))
            for i in range(len(new_weight)):
                if new_weight[i] > 0.999:
                    new_weight[i] = 0.999
                    for j in range(len(new_weight)):
                        if j != i:
                            new_weight[j] += 0.001 / (len(new_weight) - 1)
                            break
                elif new_weight[i] < 0.001:
                    new_weight[i] = 0.001
                    for j in range(len(new_weight)):
                        if j != i:
                            new_weight[j] -= 0.001 / (len(new_weight) - 1)
                            break
            new_weight = new_weight / np.sum(new_weight)
            closest_idx = np.argmin(np.linalg.norm(normalized_population_fitness - new_weight, axis=1))
            W[closest_idx] = new_weight
        self.B = Neighbor(W, self.T)
        return population_fitness, W

    def Chromo_setup(self):
        self.os_list = []
        for i in range(len(self.args.O_num)):
            self.os_list.extend([i for _ in range(self.args.O_num[i])])
        self.half_len_chromo =len(self.os_list)
        self.ms_list =[]
        self.J_site =[]
        for i in range(len(self.args.Processing_Machine)):
            for j in range(len(self.args.Processing_Machine[i])):
                self.ms_list.append(len(self.args.Processing_Machine[i][j]))
                self.J_site.append((i ,j))

    def Chromo_setup_0(self):
        self.os_list = []
        for i in range(len(self.args.O_num)):
            self.os_list.extend([i for _ in range(self.args.O_num[i])])

    def random_initial(self):
        for i in range(self.Pop_size):
            Pop_i =[]
            random.shuffle(self.os_list)
            Pop_i.extend(copy.copy(self.os_list))
            ms =[]
            for i in self.ms_list:
                ms.append(random.randint(0 , i -1))
            Pop_i.extend(ms)
            Pop_i =Popi(self.Ihang, self.JLie,self.args ,Pop_i ,self.J_site ,self.half_len_chromo)
            if self._z==[]:
                self._z =Pop_i.fitness
            else:
                for j in range(2):
                    if self._z[j ] >Pop_i.fitness[j]:
                        self._z[j ] =Pop_i.fitness[j]
            self.Pop.append(Pop_i)


    # POX:precedence preserving order-based crossover
    def POX(self ,p1, p2):
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

    def Job_Crossover(self ,p1 ,p2):
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
        return o1 ,o2

    def swap_mutation(self ,p):
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

    def MB_mutation(self ,p1):
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
        r = random.randint(1, self.half_len_chromo)
        random.shuffle(T_r)
        R = T_r[0:r]
        for i in R:
            K, K_2 = CHS1[i], CHS2[i]
            CHS1[i], CHS2[i] = K_2, K
        return CHS1, CHS2

    def Mutation_Machine(self ,CHS):
        T_r = [j for j in range(self.half_len_chromo)]
        r = random.randint(1, self.half_len_chromo)
        random.shuffle(T_r)
        R = T_r[0:r]
        for i in R:
            O_site =self.J_site[i]
            pt =self.args.Processing_Time[O_site[0]][O_site[1]]
            pt_find =pt[0]
            len_pt =len(pt ) -1
            k , m =1 ,0
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
        P1, P2 = Popi(self.Ihang, self.JLie, self.args, p1, self.J_site, self.half_len_chromo), \
            Popi(self.Ihang, self.JLie, self.args, p2, self.J_site, self.half_len_chromo)
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
        P1, P2 = Popi(self.Ihang, self.JLie, self.args, p1, self.J_site, self.half_len_chromo), \
            Popi(self.Ihang, self.JLie, self.args, p2, self.J_site, self.half_len_chromo)
        return P1, P2

    def adaW_main(self, Ihang, Jlie):

        # ----------------------------------Initialization---------------------------------
        # to obtain Populations and weight vectors
        self.Pop_size = len(self._lambda)
        self.Pop=[]
        if self.means_m > 1:  # Used to determine whether the machine is flexible
            self.random_initial()  # random Initial

        self.B = Neighbor(self._lambda, self.T)
        EP = []  # EP is used to store non-dominated solutions found during the search
        # ----------------------------------Iteration---------------------------------
        for gi in range(self.gene_size):
            self.pc = 0.8
            self.pm = 0.05
            for i in range(len(self.Pop)):
                # Randomly select two indexes k,l from B(i)
                j = random.randint(0, self.T - 1)
                k = random.randint(0, self.T - 1)
                # generate new solution from pop[j] and pop[k] by using genetic operators
                if self.means_m > 1:
                    pop1, pop2 = self.operator_Flexible(self.Pop[self.B[i][j]], self.Pop[self.B[i][k]])
                else:
                    pop1, pop2 = self.operator_NoFlexible(self.Pop[self.B[i][j]], self.Pop[self.B[i][k]])
                is1 = False
                if Tri_Dominate(pop1, pop2):
                    y = pop1
                    is1 = True
                else:
                    y = pop2
                # update of the reference point z
                for zi in range(len(self._z)):
                    if self._z[zi] > y.fitness[zi]:
                        self._z[zi] = y.fitness[zi]
                        if zi == 2:
                            self.min_1fit = y.fitness[zi]
                            if is1:
                                self.min_pop = pop1
                            else:
                                self.min_pop = pop2
                # update of Neighboring solutions
                for bi in range(len(self.B[i])):
                    Ta = Tchebycheff(self.Pop[self.B[i][bi]], self._z, self._lambda[self.B[i][bi]])
                    Tb = Tchebycheff(y, self._z, self._lambda[self.B[i][bi]])
                    if Tb < Ta:
                        self.Pop[self.B[i][j]] = y
                        break
                # Update of EP
                if EP == []:
                    EP.append(y)
                else:
                    dominateY = False
                    _remove = []  # Remove from EP all the vectors dominated by y
                    _remove_l = []
                    for ei in range(len(EP)):
                        if Tri_Dominate(y, EP[ei]):
                            _remove.append(EP[ei])
                        elif Tri_Dominate(EP[ei], y):
                            dominateY = True
                            break
                    # add y to EP if no vectors in EP dominated y
                    if not dominateY:
                        EP.append(y)

                        for j in range(len(_remove)):
                            EP.remove(_remove[j])
            if (gi+1) % 2 == 0 and gi < self.Pop_size * 0.9:
                self.weight_update(EP)
        return EP, self._z
