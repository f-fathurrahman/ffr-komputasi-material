from sympy import *
init_printing(use_unicode=True)

εA, εB = symbols("epsilon_A epsilon_B", real=True)
E = symbols("E", real=True)
t = symbols("t", real=True)
k, a = symbols("k a")

secMat = Matrix([
    [εA - E, -t*( 1 + cos(k*a) - I*sin(k*a) )],
    [-t*( 1 + cos(k*a) + I*sin(k*a) ), εB - E]
])

dict_subs = {
    εA: 1.0,
    εB: 1.1,
    t: 0.5,
    a: 1.0,
    k: 0.0
}

import numpy as np
import matplotlib.pyplot as plt

knums = np.linspace(-2*np.pi, 2*np.pi, 101)
plt.figure(figsize=[6.0, 10.0])
plt.clf()
for knum in knums:
    dict_subs[k] = knum
    secMat_num = secMat.subs(dict_subs)
    #pprint(secMat_num)
    detsecMat = secMat_num.det()
    sols = solve(detsecMat, E)
    Nroots = len(sols)
    for i in range(Nroots):
        plt.plot([knum], sols[i], marker="o", color="b", markersize=2.0)
    print("knum = ", knum, " is done")

plt.grid(True)
plt.savefig("IMG_tb1d_diatomik_v1.pdf")
