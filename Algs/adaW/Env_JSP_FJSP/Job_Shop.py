from Env_JSP_FJSP.Job import Job
from Env_JSP_FJSP.Machine import Machine


class Job_shop:
    def __init__(self, args,Ihang, JLie):
        self.n = args.n
        self.m = args.m
        self.O_num = args.O_num
        self.PM = args.Processing_Machine
        self.PT = args.Processing_Time
        self.PQ = args.Processing_Quality
        self.PS = args.Processing_Speed
        self.reset(Ihang, JLie)

    def reset(self,Ihang, JLie):
        self.C_max = 0  # makespan
        self.energy = 0
        self.single_energy = [0] * self.m
        self.energy_statistic = [0 for _ in range(self.m)]
        self.energy_oldsta = [0 for _ in range(self.m)]
        self.load = 0  # Total load of machines
        self.single_load = [0] * self.m
        self.max_EndM = None  # the last end machine
        self.Jobs = []
        for i in range(self.n):
            Ji = Job(i, self.PM[i], self.PT[i], self.PS, self.PQ,Ihang, JLie)
            self.Jobs.append(Ji)
        self.Machines = []
        for j in range(self.m):
            Mi = Machine(j,Ihang, JLie)
            self.Machines.append(Mi)

    # decode of chs[i]
    def decode(self, job, machine):
        Ji = self.Jobs[job]  # state of current job
        # obtain processing time/start time/processing machine of current operation
        p_q, p_s = Ji.get_processing_speed()
        o_pt, last_endt, M_idx = Ji.get_next_info(machine)
        Mi = self.Machines[M_idx - 1]  # state of current machine
        start = Mi.find_start(last_endt, o_pt)  # obtatin real start time on machine
        end = start + o_pt
        self.load += o_pt
        self.single_load[M_idx-1] += o_pt
        Mi.update(start, end, p_q, [Ji.idx, Ji.cur_op])  # update machine state  [Ji.idx, Ji.cur_op]：第Ji.idx个零件（有0），第Ji.cur_op个工序（有0）
        Ji.update(start, end, Mi.idx, p_q, p_s)  # update Job state
        self.energy_oldsta[Mi.idx] = self.energy_old(Mi)
        # self.energy_statistic[Mi.idx] = self.energy_rule(Mi,Ji)
        self.energy = sum(self.energy_oldsta)

        if end > self.C_max:  # update makespan
            self.C_max = end
            self.max_EndM = Mi
    def energy_old(self,mi):
        total_energy = 0
        if len(mi.start) == 1:
            for i in range(len(mi.start)):
                process_t = mi.end[i] - mi.start[i]
                p_n = mi._on[i]
                u_e = self.Jobs[p_n[0]].o_e[p_n[1]]
                total_energy += process_t * u_e
        elif len(mi.start) > 1:
            for i in range(len(mi.start)):
                process_t = mi.end[i] - mi.start[i]
                p_n = mi._on[i]
                u_e = self.Jobs[p_n[0]].o_e[p_n[1]]
                total_energy += process_t * u_e
                if i < len(mi.start)-1:
                    window_t = mi.start[i+1] - mi.end[i]
                    w_e = mi.standbye
                    total_energy += window_t * w_e
        total_energy = total_energy + mi.updown
        return total_energy
