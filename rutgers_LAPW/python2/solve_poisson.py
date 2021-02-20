import numpy as np
from math import pi
import weave

def solve_poisson(Zq, R, rho):
    """Given the input density rho, calculates the Hartree potential
    The boundary conditions used are U(0)=0 and U(S)=Zq. The boundary condition at S is only a constant shift
    of potential which is later readjusted by choosing MT zero. So, this boundary condition is not very relevant.
    """
    codeNumerovInh="""
        double h2 = dx*dx;
        double h12 = h2/12;
        
        double w0 = Solution(0)-h12*U(0);
        double Ux = U(1);
        double w1 = Solution(1)-h12*Ux;
        double Phi = Solution(1);
        
        double w2;
        for (int i=2; i<Nmax; i++){
          w2 = 2*w1 - w0 + h2*Ux;
          w0 = w1;
          w1 = w2;
          Ux = U(i);
          Phi = w2+h12*Ux;
          Solution(i) = Phi;
        }
     """
    U = np.array([-4*pi*r*rho[i] for i,r in enumerate(R)])
    Nmax = len(R)
    dx = float( (R[-1]-R[0])/(len(R)-1.) )
    Solution = np.zeros(len(R), dtype=float)
    
    Solution[0]=0
    Solution[1]=(R[1]-R[0]) # Boundary condition for U_H=V_H/r
    
    weave.inline(codeNumerovInh, ['U', 'Nmax', 'dx', 'Solution'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    # adding homogeneous solution to satisfay boundary conditions: U(0)=0, U(infinity)=Z
    alpha = (Zq - Solution[-1])/R[-1]
    Solution += alpha*R
    return Solution
