from sympy import *
init_printing(use_unicode=True)

k = symbols("k", real=True)
a = symbols("a", real=True, positive=True)
m1, m2 = symbols("m1 m2", real=True, positive=True)
κ1, κ2 = symbols("kappa1 kappa2", real=True, positive=True)
ω = symbols("omega", real=True)

M = Matrix([
    [2*((κ1 + κ2) - κ2*cos(k*a))/m1, -κ1*(1 + cos(k*a) - I*sin(k*a))],
    [-κ1*(1 + cos(k*a) + I*sin(k*a)), 2*((κ1 + κ2) - κ2*cos(k*a)) - m2*ω**2]
])

eq1 = M.det()
dict_num = {
    m1: 1.0,
    m2: 1.0,
    κ1: 1.0,
    κ2: 1.0,
    k: 0.0,
    a: 1.0
}


