from sympy import symbols,Rational,solve,simplify,N

alpha,Z = symbols('alpha Z')
atoms = ['H','He','Li','Be']
exacts = [-0.5,-2.9037,-7.47806,-14.66736]

E_kin = Rational(1,2)*alpha**2
E_en = -Z*alpha
E_ee = Rational(5,8)*alpha

def E(Z,alpha):
    Nel = Z
    return Z*(E_kin-Z*alpha)+Rational(Nel*(Nel-1),2)*E_ee

print('Theoretical bosonic energies for 1S Ansatz, nucleus mass is infinite')
for i,atom in enumerate(atoms):
    Z = i+1
    ex = exacts[i]
    EE = simplify(E(Z,alpha))
    dE = EE.diff(alpha)
    alpha_min = solve(dE)[0]
    Emin =  EE.subs(alpha,alpha_min)
    print("{:1} {:>2} {:>30} minimum at alpha: {:<6} = {:<8.5}  E min: {:10}={:10.5f} Exact E:{:10.5f}"\
          .format(Z,atom,str(EE),str(alpha_min),N(alpha_min),str(Emin),N(Emin),ex))

