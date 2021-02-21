from sympy import *
init_printing(use_unicode=True)

k = symbols("k", real=True)
a = symbols("a", real=True, positive=True)
m1, m2 = symbols("m1 m2", real=True, positive=True)
κ1, κ2 = symbols("kappa1 kappa2", real=True, positive=True)
ω = symbols("omega", real=True)

# Using the secular equation given in the text
eq1 = m1*m2*ω**4 - (m1 + m2)*(2*κ1 + 4*κ2*sin(k*a/2)**2)*ω**2 + \
(4*κ1**2 + 16*κ1*κ2)*sin(k*a/2)**2 + 16*κ2**2*sin(k*a/2)**4
dict_num = {
    m1: 1.0,
    m2: 1.1,
    κ1: 2.0,
    κ2: 0.1,
    k: 0.0,
    a: 1.0
}


import numpy as np
import matplotlib.pyplot as plt

NkPlot = 100
k_nums = np.linspace(-2*np.pi, 2*np.pi, NkPlot)
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
        if abs(im(x)) > 1e-10:
            print("WARNING: Imaginary part detected")
    sols = np.array(solsn)
    sols = sols[sols >= 0]
    Nelems = len(sols)
    for i in range(Nelems):
        plt.plot([k_n], sols[i], marker="o", color="b", markersize=0.5)
    print("done k_n = ", k_n)

plt.grid(True)
plt.savefig("IMG_Simon_10_4_v2.pdf")