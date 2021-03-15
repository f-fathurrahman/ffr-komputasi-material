from sympy import *
init_printing(use_unicode=True)

k = symbols("k", real=True)
a = symbols("a", real=True, positive=True)
m1, m2 = symbols("m1 m2", real=True, positive=True)
κ1, κ2 = symbols("kappa1 kappa2", real=True, positive=True)
ω = symbols("omega", real=True)


dynMat = Matrix([
    [(κ1 + κ2)/m1, (-κ2 - κ1*(cos(k*a) + I*sin(k*a)))/m1 ],
    [(-κ2 - κ1*(cos(k*a) - I*sin(k*a)))/m2 , (κ1 + κ2)/m2]
])

dict_num = {
    m1: 1.0,
    m2: 1.0,
    κ1: 1.0,
    κ2: 2.1,
    k: 0.0,
    a: 1.0
}
dynMat_num = dynMat.subs(dict_num)
pprint(dynMat_num)


import numpy as np
import matplotlib.pyplot as plt

NkPlot = 101
k_nums = np.linspace(-np.pi, np.pi, NkPlot)
plt.clf()
for k_n in k_nums:
    dict_num[k] = k_n
    dynMat_num = dynMat.subs(dict_num)
    res = dynMat_num.eigenvects()
    solsn = []
    for x in res:
        solsn.append(re(x[0]))
    sols = np.array(solsn)
    sols = sols[sols >= 0]
    Nelems = len(sols)
    for i in range(Nelems):
        plt.plot([k_n], sqrt(sols[i]), marker="o", color="b", markersize=0.5)
    print("done k_n = ", k_n)

plt.grid(True)
plt.savefig("IMG_fonon_eigen.pdf")
