import json
import configparser

class Machine:
    def __init__(self, idx,Ihang, JLie):
        self.idx = idx
        self.start = []
        self.end = []
        self.p_ql = []
        self.init_energy(idx,Ihang, JLie)
        self._on = []
        self.ced = []
        self.calenergy_ok = []
        self.old_start = []
        self.old_end = []

    def init_energy(self,idx,Ihang, JLie):
        if Ihang < 9:
            f = 'Mk0' + str(Ihang + 1) + '.pkl'
        else:
            f = 'Mk' + str(Ihang + 1) + '.pkl'
        job_dic, ma_dic = self.get_inputInfo(f, Ihang, JLie)
        self.safet = ma_dic["safet"][idx]
        self.standbye = ma_dic["standbye"][idx]
        self.updown = ma_dic["updown"][idx]
        self.main_safe = ma_dic["main_safe"][idx]
        self.energyt = ma_dic["energyt"][idx]



    def update(self, s, e, p_q, Job):
        self.start.append(s)
        self.start.sort()
        self.end.append(e)
        self.end.sort()
        self.old_start.append(s)
        self.old_start.sort()
        self.old_end.append(e)
        self.old_end.sort()
        idx = self.start.index(s)
        self._on.insert(idx, Job)
        self.p_ql.append(p_q)
        self.ced.append(True)
        self.calenergy_ok.append(0)


    def find_start(self, s, o_pt):
        if self.end == []:
            return max(s, 0)
        else:
            if s > self.end[-1]:
                return s
            else:
                o_s = self.end[-1]
                l = len(self.end) - 2
                while l >= 0:
                    if s + o_pt > self.start[l + 1]:
                        break
                    if self.end[l] > s and self.end[l] + o_pt <= self.start[l + 1] and (self.calenergy_ok[l] == 0 and self.calenergy_ok[l+1] == 0):
                        o_s = self.end[l]
                    elif self.end[l] < s and s + o_pt <= self.start[l + 1] and (self.calenergy_ok[l] == 0 and self.calenergy_ok[l+1] == 0):
                        o_s = s
                    l -= 1
                return o_s

    def re_update(self, s, e, Job):
        self.start.append(s)
        self.start.sort()
        self.end.append(e)
        self.end.sort()
        idx = self.start.index(s)
        self._on.insert(idx, Job)

    def get_inputInfo(self,key,Ihang, JLie):
        config = configparser.ConfigParser()
        config.read(f'E:\PaperCode\Result\config\Brandimarte\parameter{Ihang}.ini')
        ma_info = config.get(key, 'ma_dic')
        ma_dic = json.loads(ma_info)
        job_info = config.get(key, 'job_dic')
        job_dic = json.loads(job_info)
        return job_dic,ma_dic
