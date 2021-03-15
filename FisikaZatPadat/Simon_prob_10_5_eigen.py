from sympy import *
init_printing(use_unicode=True)

k = symbols("k", real=True)
a = symbols("a", real=True, positive=True)
m1, m2, m3 = symbols("m1 m2 m3", real=True, positive=True)
κ1, κ2, κ3 = symbols("kappa1 kappa2 kappa3", real=True, positive=True)
ω = symbols("omega", real=True)

dynMat = Matrix([
    [(κ3 + κ1)/m1, κ1/m1, κ3*( cos(k*a) - I*sin(k*a) )/m1],
    [κ1/m2,  (κ1 + κ2)/m2, κ2/m2],
    [κ3*( cos(k*a) + I*sin(k*a) )/m3, κ2/m3, (κ2 + κ3)/m3]
])

dict_num = {
    m1: 1.0,
    m2: 1.1,
    m3: 1.2,
    κ1: 1.0,
    κ2: 1.1,
    κ3: 1.5,
    k: 1.0,
    a: 1.0
}
