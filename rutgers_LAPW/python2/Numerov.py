import numpy as np
import weave

########################################################
# Routines for solving the ODE problem and Poisson EQ. #
########################################################
def Numerov(F, dx, f0=0.0, f1=1e-3):
    codeNumerov="""
      double h2 = dx*dx;
      double h12 = h2/12;
      
      double w0 = (1-h12*F(0))*Solution(0);
      double Fx = F(1);
      double w1 = (1-h12*Fx)*Solution(1);
      double Phi = Solution(1);
      
      double w2;
      for (int i=2; i<Nmax; i++){
        w2 = 2*w1 - w0 + h2*Phi*Fx;
        w0 = w1;
        w1 = w2;
        Fx = F(i);
        Phi = w2/(1-h12*Fx);
        Solution(i) = Phi;
      }
    """
    Nmax = len(F)
    dx = float(dx)
    Solution = np.zeros(Nmax, dtype=float)
    Solution[0] = f0
    Solution[1] = f1
    weave.inline(codeNumerov, ['F', 'Nmax', 'dx', 'Solution'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return Solution

def NumerovGen(F, U, dx, f0=0.0, f1=1e-3):
    codeNumerov="""
      double h2 = dx*dx;
      double h12 = h2/12;
      
      double w0 = Solution(0)*(1-h12*F(0))-h12*U(0);
      double Fx = F(1);
      double Ux = U(1);
      double w1 = Solution(1)*(1-h12*Fx)-h12*Ux;
      double Phi = Solution(1);
      
      double w2;
      for (int i=2; i<Nmax; i++){
        w2 = 2*w1 - w0 + h2*(Phi*Fx+Ux);
        w0 = w1;
        w1 = w2;
        Fx = F(i);
        Ux = U(i);
        Phi = (w2+h12*Ux)/(1-h12*Fx);
        Solution(i) = Phi;
      }
    """
    Nmax = len(F)
    dx = float(dx)
    Solution = np.zeros(Nmax, dtype=float)
    Solution[0] = f0
    Solution[1] = f1
    weave.inline(codeNumerov, ['F', 'U', 'Nmax', 'dx', 'Solution'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return Solution
