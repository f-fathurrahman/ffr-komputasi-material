import sys
import numpy as np

def read_abinit_eig(filename):
    
    f = open(filename)
    
    line = f.readline()
    strs = line.split()
    
    print(strs)
    unit = strs[2]
    Nkpt = int(strs[6])
    print("Nkpt = ", Nkpt)

    line = f.readline()
    line = line.replace(",","")
    strs = line.split()

    print(strs)
    Nband = int(strs[3])
    print("Nband = ", Nband)

    ebands = np.zeros((Nband,Nkpt))

    for ik in range(Nkpt):
        line = f.readline()
        strs = line.split()
        line = f.readline()



read_abinit_eig(sys.argv[1])