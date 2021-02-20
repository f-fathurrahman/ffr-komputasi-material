import numpy as np
import weave
from math import pi, exp

def calc_MT_density(mu, Ek, wkp, w0, w1, w2, Psi_l, Psip_l, beta=50.):
    """Given the coefficients Eqs.58-61 on page 30, it computes the valence charge
       given the chemical potential mu. The Eq. 36 on page 31 descibes the algorithm.
       Eq.36 shows that the radial solution of the SCH equation is necessary together with its energy derivative
    """
    code_sum="""
    double dw0, dw1, dw2; dw0=dw1=dw2=0;
    for (int p=0; p<ek.shape()[0]; p++){
        double x = beta*(ek(p)-mu);
        double ferm;
        ferm = abs(x) < 100 ? 1/(exp(x)+1) : (x<0 ? 1 : 0);
        dw0 += w0k(p)*ferm;
        dw1 += w1k(p)*ferm;
        dw2 += w2k(p)*ferm;
    }
    py::tuple results(3);
    results[0] = dw0; results[1] = dw1; results[2] = dw2;
    return_val = results;
    """
    nlmax = len(w0[0])
    wgh = np.zeros((nlmax,3), dtype=float)
    for l in range(nlmax):
        for ik in range(np.shape(Ek)[0]):
            w0k = np.array(w0[ik][l,:])
            w1k = np.array(w1[ik][l,:])
            w2k = np.array(w2[ik][l,:])
            ek = np.array(Ek[ik])
            dws = weave.inline(code_sum, ['mu', 'ek', 'beta', 'w0k', 'w1k', 'w2k'],type_converters=weave.converters.blitz, compiler = 'gcc')
            wgh[l,:] += np.array(dws) * wkp[ik]
        #print "%3d %20.10f %20.10f %20.10f" % (l, wgh[l,0], wgh[l,1], wgh[l,2])

    # Implements Eq.63 on page 31.
    nR = len(Psi_l[0])
    MTRho = np.zeros(nR, dtype=float)
    for l in range(nlmax):
        for ir in range(nR):
            MTRho[ir] += wgh[l,0]*Psi_l[l][ir]**2 + wgh[l,1]*Psi_l[l][ir]*Psip_l[l][ir] + wgh[l,2]*Psip_l[l][ir]**2
            
    MTRho *= 2/(4*pi)  # 2 due to spin

    return MTRho