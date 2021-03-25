import weave
import numpy as np

def CRHS(E, l, R, Veff):
    "RHS for solving the Schroedinger equations by Numerov. To achive sufficient speed, uses C++."
    codeRHS="""
        for (int i=0; i<N; i++){
           RHS(i) = 2*( -E + 0.5*l*(l+1)/(R(i)*R(i)) + Veff(i) );
        }
    """
    N = len(R)
    RHS = np.zeros(len(R), dtype=float)
    weave.inline(codeRHS, ['N', 'E', 'l', 'R', 'Veff', 'RHS'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return RHS
