from Algorithms.Params import get_args
import datetime
from Env_JSP_FJSP.Speed_Select import Speed_T_Matrix, Processing_Quality
from Algorithms.Algorithm import *
from adaW.Algorithms.utils import *


def Algo_Solver(f, Ihang, JLie):
    file = r'E:\PaperCode\data\MK' + '/' + f#Brandimarte
    n, m, PT, MT, ni = Instance(file)
    PS = Speed_T_Matrix
    PQ = Processing_Quality
    mm = 2
    args = get_args(n, m, PT, MT, PS, PQ, ni, mm)
    Algo = Algorithms(args, Ihang, JLie)
    adaWEP,M_z,EP_out_all = Algo.mogwo_main()
    input_json(Ihang, JLie, adaWEP)


if __name__=="__main__":
    now_start = datetime.datetime.now()
    print("-----------------------------------------------")
    print("start：", now_start)
    for k in range(10):
        print('Time: ', k + 1)
        for j in range(10):
            print('time: ', j+1)
            for i in range(10):
                if i<9:f = 'Mk0'+str(i+1)+'.pkl'
                else:f = 'Mk'+str(i+1)+'.pkl'
                print(f)
                Algo_Solver(f,i,j)
    now_end = datetime.datetime.now()
    print("-----------------------------------------------")
    print("end：", now_end)
