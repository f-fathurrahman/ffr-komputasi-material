from sympy import *
init_printing(use_unicode=True)

k = symbols("k", real=True)
a = symbols("a", real=True, positive=True)
m1, m2 = symbols("m1 m2", real=True, positive=True)
κ1, κ2 = symbols("kappa1 kappa2", real=True, positive=True)
ω = symbols("omega", real=True)

secMat = Matrix([
    [(κ1 + κ2) - m1*ω**2, -κ2 - κ1*(cos(k*a) + I*sin(k*a)) ],
    [-κ2 - κ1*(cos(k*a) - I*sin(k*a)), (κ1 + κ2) - m2*ω**2]
])

eq1 = secMat.det()
dict_num = {
    m1: 1.0,
    m2: 1.0,
    κ1: 1.0,
    κ2: 1.5,
    k: 0.0,
    a: 1.0
}


import numpy as np
import matplotlib.pyplot as plt

NkPlot = 101
k_nums = np.linspace(-np.pi, np.pi, NkPlot)
plt.clf()
for k_n in k_nums:
    dict_num[k] = k_n
    print(dict_num)
    eq1num = eq1.subs(dict_num)
    print(eq1num)
    sols = roots(eq1num) # use roots instead of solve
    solsn = []
    for x in sols.keys():
        solsn.append(re(x))
    sols = np.array(solsn)
    sols = sols[sols >= 0]
    Nelems = len(sols)
    for i in range(Nelems):
        plt.plot([k_n], sols[i], marker="o", color="b", markersize=0.5)
    print("done k_n = ", k_n)

plt.grid(True)
plt.xlabel("k")
plt.ylabel("ω")
plt.savefig("IMG_fonon_diatomik_v1.pdf")
