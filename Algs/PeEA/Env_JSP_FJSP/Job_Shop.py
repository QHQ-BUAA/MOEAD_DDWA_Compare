from Env_JSP_FJSP.Job import Job
from Env_JSP_FJSP.Machine import Machine

class Job_shop:
    def __init__(self,args,job_dic, ma_dic):
        self.n= args.n
        self.m=args.m
        self.O_num=args.O_num
        self.PM = args.Processing_Machine
        self.PT = args.Processing_Time
        self.reset(job_dic, ma_dic)

    def reset(self,job_dic, ma_dic):
        self.C_max = 0      #makespan
        self.energy = 0
        self.max_EndM=None  # the last end machine
        self.energy_oldsta = [0 for _ in range(self.m)]
        self.single_load= [0] * self.m    # load of each machine
        self.Jobs=[]
        for i in range(self.n):
            Ji=Job(i,self.PM[i],self.PT[i],job_dic)
            self.Jobs.append(Ji)
        self.Machines=[]
        for j in range(self.m):
            Mi=Machine(j, ma_dic)
            self.Machines.append(Mi)

    # decode of chs[i]
    def decode(self,Job,Machine):
        Ji=self.Jobs[Job]  #state of current job
        # obtain processing time/start time/processing machine of current operation
        o_pt, s,M_idx = Ji.get_next_info(Machine)
        Mi=self.Machines[M_idx-1]   #state of current machine
        start=Mi.find_start(s,o_pt)     # obtatin real start time on machine
        end=start+o_pt
        self.single_load[Mi.idx]+=o_pt
        Mi.update(start, end, [Ji.idx, Ji.cur_op])
        Ji.update(start,end,Mi.idx)     #update Job state
        self.energy_oldsta[Mi.idx] = self.energy_old(Mi)
        self.energy = sum(self.energy_oldsta)
        if end>self.C_max:  # update makespan
            self.C_max=end
            self.max_EndM=Mi
        self.max_load = max(self.single_load)  #update max_load of machine

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

    # endregion