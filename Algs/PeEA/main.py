from Algorithms.Algorithm import *
from Algorithms.Params import get_args
import datetime
from PeEA.Algorithms.utils import *


def Algo_Solver(Ihang,Jlie,f):
    file = r'E:\PaperCode\data\Brandimarte' + '/' + f
    n, m, PT, MT, ni = Instance(file)
    mm = 2
    job_dic,ma_dic = get_inputInfo(f,Ihang,Jlie)
    args = get_args(n, m, PT, MT, ni, mm)
    Algo = Algorithms(args,job_dic,ma_dic)
    PeEA_EP = Algo.PeEA_main(Ihang, Jlie)
    input_json(Ihang, Jlie, PeEA_EP)


if __name__=="__main__":
    now_start = datetime.datetime.now()
    print("-----------------------------------------------")
    print("start：", now_start)
    for j in range(10):
        print('time: ', j+1)
        for i in range(10):
            if i<9:f = 'Mk0'+str(i+1)+'.pkl'
            else:f = 'Mk'+str(i+1)+'.pkl'
            print('instance: ', f)
            Algo_Solver(i, j, f)
    now_end = datetime.datetime.now()
    print("-----------------------------------------------")
    print("end：", now_end)
