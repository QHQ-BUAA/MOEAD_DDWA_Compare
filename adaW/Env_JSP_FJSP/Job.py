import random
from random import choice
from Env_JSP_FJSP.Speed_Select import Operation_Energy
import json
import configparser


class Job:
    def __init__(self, idx, processing_machine, processing_time, processing_speed, processing_quality,Ihang, JLie):
        self.idx = idx
        self.p_m = processing_machine
        self.p_t = processing_time
        self.p_q = processing_quality
        self.p_s = processing_speed
        self.O_e = Operation_Energy
        self.p_qlist = []
        self.p_slist = []
        self.end = []
        self.cur_op = 0
        self.cur_pt = None
        self._on = []
        self.start = []
        self.endt = 0
        self.old_oe = []

        if Ihang < 9:
            f = 'Mk0' + str(Ihang + 1) + '.pkl'
        else:
            f = 'Mk' + str(Ihang + 1) + '.pkl'
        job_dic, ma_dic = self.get_inputInfo(f, Ihang, JLie)
        self.p_qlist = job_dic['pq_list'][idx]
        self.p_slist = job_dic['ps_list'][idx]
        self.o_e = job_dic['job_oe'][idx]
        for i in self.o_e:
            self.old_oe.append(i)

    def get_next_info(self, Machine):
        m_idx = self.p_m[self.cur_op][Machine]
        self.cur_pt = self.p_t[self.cur_op][Machine]
        return self.cur_pt, self.endt, m_idx

    def update(self, s, e, m_idx, p_q, p_s):
        self.cur_op += 1
        self.start.append(s)
        self.end.append(e)
        self._on.append(m_idx)




    def get_processing_speed(self):
        random_index = random.randint(0, len(self.p_q) - 1)
        processing_quailty = self.p_q[random_index]
        if processing_quailty > 6:
            processing_speed_list = self.p_s[5:]
            processing_speed = choice(processing_speed_list)
        else:
            processing_speed_list = self.p_s[random_index:]
            processing_speed = choice(processing_speed_list)
        return processing_quailty, processing_speed

    def get_inputInfo(self,key,Ihang, JLie):
        config = configparser.ConfigParser()
        config.read(f'E:\PaperCode\Result\config\Brandimarte\parameter{Ihang}.ini')
        ma_info = config.get(key, 'ma_dic')
        ma_dic = json.loads(ma_info)
        job_info = config.get(key, 'job_dic')
        job_dic = json.loads(job_info)
        return job_dic,ma_dic
