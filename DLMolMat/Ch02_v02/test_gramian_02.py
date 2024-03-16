from sympy import *

ε = symbols("epsilon")
δ = symbols("delta")
X = Matrix([
    [1, 1 + ε, 1 + 2*ε],
    [1, 1 + δ, 1 + 2*δ]
])
X[0,0] += 3

G = X.T * X
pprint(G)
print("Determinant of G = ")
pprint(G.det())

