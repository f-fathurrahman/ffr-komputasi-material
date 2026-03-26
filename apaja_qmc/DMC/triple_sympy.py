from sympy import *
from sympy.physics.quantum.operator import *

tau = Symbol('tau')
A = Operator('A')
B = Operator('B')
C = Operator('C')

print('Exact formula expanded to 4th order:')
print('='*40)
H = A+B+C
serH = series(exp(-tau*H),tau,0,3)
print('exp(-tau*H) = exp(-tau*(A+B+C)) = ')
pprint(serH)

print('2nd order formulas:')
print('='*40)
expABCBA = exp(-tau/2*A)*exp(-tau/2*B)*exp(-tau*C)*exp(-tau/2*B)*exp(-tau/2*A)
print(str(expABCBA)+' = ')
serABCBA =  series(expABCBA,tau,0,3)
pprint(serABCBA)


print("Warning: SymPy expands all exp(-tau*A) to series without checking convergence")
