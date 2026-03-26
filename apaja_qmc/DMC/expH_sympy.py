from sympy import *
from sympy.physics.quantum.operator import *

tau = Symbol('tau')
#H = HermitianOperator('H')
T = Operator('T')
V = Operator('V')

print('Exact formula expanded to 4th order:')
print('='*40)
H = T+V
serH = series(exp(-tau*H),tau,0,4)
print('exp(-tau*H) = ')
pprint(serH)

print('2nd order formulas:')
print('='*40)
expTVT = exp(-tau/2*T)*exp(-tau*V)*exp(-tau/2*T)
serTVT = series(expTVT,tau,0,4)
print('exp(-tau/2*T)*exp(-tau*V)*exp(-tau/2*T) = ')
pprint(serTVT)


expVTV = exp(-tau/2*V)*exp(-tau*T)*exp(-tau/2*V)
serVTV = series(expVTV,tau,0,4)
print('exp(-tau/2*V)*exp(-tau*T)*exp(-tau/2*V) = ')
pprint(serVTV)

print("Warning: SymPy expands exp(-tau*T) to series without checking it's convergence")
