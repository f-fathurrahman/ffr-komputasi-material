from sympy import *
init_printing(use_unicode=True)

k = symbols("k", real=True)
a = symbols("a", real=True, positive=True)
m1, m2, m3 = symbols("m1 m2 m3", real=True, positive=True)
κ1, κ2, κ3 = symbols("kappa1 kappa2 kappa3", real=True, positive=True)
ω = symbols("omega", real=True)

M = Matrix([
    [κ3 + κ1 - m1*ω**2,            κ1,                κ3*( cos(k*a) - I*sin(k*a) )],
    [κ1,                           κ1 + κ2 - m2*ω**2, κ2],
    [κ3*( cos(k*a) + I*sin(k*a) ), κ2,                κ2 + κ3 - m3*ω**2]
])

eq1 = M.det()
dict_num = {
    m1: 1.0,
    m2: 1.1,
    m3: 1.1,
    κ1: 1.0,
    κ2: 1.1,
    κ3: 1.5,
    k: 0.0,
    a: 1.0
}

import numpy as np
import matplotlib.pyplot as plt

NkPlot = 100
k_nums = np.linspace(-2*np.pi, np.pi, NkPlot)
plt.clf()
for k_n in k_nums:
    dict_num[k] = k_n
    print(dict_num)
    eq1num = eq1.subs(dict_num)
    print(eq1num)
    sols = np.array( list(roots(eq1num).keys()) ) # use roots instead of solve
    sols = sols[sols >= 0]
    Nelems = len(sols)
    for i in range(Nelems):
        plt.plot([k_n], sols[i], marker="o", color="b", markersize=0.5)
    print("done k_n = ", k_n)

plt.grid(True)
plt.savefig("IMG_Simon_10_5_v1.pdf")