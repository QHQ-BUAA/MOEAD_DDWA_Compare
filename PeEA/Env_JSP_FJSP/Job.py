
class Job:
    def __init__(self, idx, processing_machine, processing_time,job_dic):
        self.idx = idx
        self.processing_machine = processing_machine
        self.processing_time = processing_time
        self.end = []
        self.cur_op = 0
        self.cur_pt = None
        self._on = []
        self.start = []
        self.endt = 0
        self.pq_list = job_dic['pq_list'][idx]
        self.ps_list = job_dic['ps_list'][idx]
        self.o_e = job_dic['job_oe'][idx]
        self.old_oe = []
        for oe in self.o_e:
            self.old_oe.append(oe)

    def get_next_info(self, Machine):
        m_idx = self.processing_machine[self.cur_op][Machine]
        self.cur_pt = self.processing_time[self.cur_op][Machine]
        return self.processing_time[self.cur_op][Machine], self.endt, m_idx

    def update(self, s, e, m_idx):
        self.endt = e
        self.cur_op += 1
        self.start.append(s)
        self.end.append(e)
        self._on.append(m_idx)