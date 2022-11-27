import time
from vasprun import Vasprun
from vasp2igor import Vasp2Igor


if __name__ == "__main__":
    start_time = time.time()
    v = Vasp2Igor('ignore/Cu.xml') 
    # print((2*3.14)/v.a)
    # # print(v.N1)
    # # print(v.N2)
    # # print(v.N3) 
    v.writefile("main.BND") 

    print("--- %s seconds ---" %(time.time() - start_time))