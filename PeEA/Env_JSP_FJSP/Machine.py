import pandas as pd

class Machine:
    def __init__(self, idx, ma_dic):
        self.idx = idx
        self.start = []
        self.end = []
        self._on = []
        self.old_start = []
        self.old_end = []
        self.init_energy(ma_dic,idx)

    def init_energy(self,ma_dic,idx):
        self.safet = ma_dic['safet'][idx]
        self.standbye = ma_dic['standbye'][idx]
        self.updown = ma_dic['updown'][idx]
        self.main_safe = ma_dic['main_safe'][idx]
        self.energyt = ma_dic['energyt'][idx]


    def update(self, s, e, Job):
        self.old_start.append(s)
        self.old_start.sort()
        self.old_end.append(e)
        self.old_end.sort()
        self.start.append(s)
        self.start.sort()
        self.end.append(e)
        self.end.sort()
        idx = self.start.index(s)
        self._on.insert(idx, Job)

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
                    if self.end[l] > s and self.end[l] + o_pt <= self.start[l + 1]:
                        o_s = self.end[l]
                    elif self.end[l] < s and s + o_pt <= self.start[l + 1]:
                        o_s = s
                    l -= 1
                return o_s