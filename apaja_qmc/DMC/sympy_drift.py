#
# Drift-diffusion to 2nd order in tau
# 
from sympy import *

tau, sqrt_tau, D = symbols('tau sqrt_tau D', positive=True)
x = symbols('x')
eta1, eta2 = symbols('eta1 eta2')
F = Function('F')

# RK2
def xp_drift(tau):   
    #y = x + D*tau*F(x)
    #return x + D*tau/2*(F(x)+F(y))
    y = x + D*tau/2*F(x)
    return x + D*tau*F(y)

xp = series(xp_drift(tau), tau, 0, n=3)

print()
print('RK2 DRIFT')
print(20*'=')

print("x' = ")
pprint(xp.doit())  # doit needed to simplify derivative 

vel = diff(xp,tau).doit() # doit needed to simplify derivative 
pprint(" dx'/dtau = ")
pprint(vel)

# DF(x')
DF = series(D*F(xp), tau, 0, 2).simplify()

pprint(" DF(x') = ")
pprint(DF)

print('Comparison:')
if vel==DF:
    print("dx'/dtau = DF(x') to 2nd order in tau") 
else:
    print("dx'/dtau is not  DF(x') to 2nd order in tau")
