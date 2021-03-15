from sympy import *
init_printing(use_unicode=True)

εA, εB = symbols("epsilon_A epsilon_B", real=True)
E = symbols("E", real=True)
t = symbols("t", real=True)
k, a = symbols("k a")

mat = Matrix([
    [εA, -t*( 1 + cos(k*a) - I*sin(k*a) )],
    [-t*( 1 + cos(k*a) + I*sin(k*a) ), εB]
])

res = mat.eigenvects()

# res[0][2] => 1st eigenvector