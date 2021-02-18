import numpy as np
from math import sqrt, fabs, pi
from scipy import special, misc
import weave
import scipy

def solve_eigensystem(k, Km, Olap_I, Enu, logDer, RMuffinTin, Vol, VKSi=0):
    """The main part of LAPW algorithm: Implements valence H[K,K'] and O[K,K'] and diagonalizes them.
       Implements all equations on page 26 and page 30.
       The output are energy bands, eigenvectors and weight functions which can be used to compute
       electronic charge in real space.
    """
    def dlog_bessel_j(lmax, x):
        """Calculates logarithmic derivative of the spherical bessel functions
           It returns three quantities:
             (x*d/dx log(j_l(x)),  j_l(x), the product of the first two)
           for l up to lmax
           The last entry is important for singular cases: when x is zero of bessel function. In this case
           the logarithmic derivative is diverging while j*dlog(j(x))/dx is not
        """
        if (fabs(x)<1e-5):
            return [(l, x**l/misc.factorial2(2*l+1), l*x**l/misc.factorial2(2*l+1)) for l in range(lmax+1)]
        else:
            (jls, djls) = special.sph_jn(lmax,x) # all jl's and derivatives for l=[0,...lmax]
            return [(x*djls[l]/jls[l], jls[l], x*djls[l]) for l in range(lmax+1)]


    # Here we prepare some quantities which depend only on one reciprocal vector K
    # These quantities are defined on pg.26 are are
    #   omegal = \omega_{l,K}
    #   C1    -- part of C^{(2)}_l
    #   PP    -- <\cdot{psi}|\cdot{psi}>
    omegal = np.zeros((len(Km), len(Enu)), dtype=float)
    C1 = np.zeros((len(Km), len(Enu)), dtype=float)
    PP = np.array([logDer[l][4] for l in range(len(Enu))])
    for iK,K in enumerate(Km):
        Dl_jl = dlog_bessel_j(len(Enu)-1, sqrt(np.dot(k+K,k+K))*RMuffinTin)
        for l in range(len(Enu)):
            (Psi, Psip, dlogPsi, dlogPsip, PsipPsip) = logDer[l]
            (Dl, jl, jlDl) = Dl_jl[l]
            omegal[iK,l] = -Psi/Psip*(Dl-dlogPsi)/(Dl-dlogPsip)
            #C1[iK,l] = sqrt(4*pi*(2*l+1)/Vol)*jl/(Psi+omegal[iK,l]*Psip) # This is less stable, but equivalent
            #  The above formula is singular for zeros of spherical bessel function. This can happen in Gamma point due to symmetry.
            C1[iK,l] = sqrt(4*pi*(2*l+1)/Vol)*(jlDl-jl*dlogPsip)/(Psi*(dlogPsi-dlogPsip))
            
    # This part of code is too slow in Python, hence was recoded in C++
    # It computes the argument of the Legendre polynomial (see Eq.48, pg.26)
    #     argums(K,K') = (k+K)*(k+K')/(|k+K| * |k+K'|)
    #     qv(K)  = k+K
    #     qvs(K) = |k+K|
    # where K is reciprocal vector and k is momentum vector.
    codeArg="""
         for (int iK=0; iK<Km.shape()[0]; iK++){
            for (int i=0; i<3; i++) qv(iK,i) = Km(iK,i)+k(i);
            double dsum=0;
            for (int i=0; i<3; i++) dsum += qv(iK,i)*qv(iK,i);
            qvs(iK) = sqrt(fabs(dsum));
         }
         for (int iK=0; iK<Km.shape()[0]; iK++){
            for (int jK=0; jK<Km.shape()[0]; jK++){
               double qvqv=0;
               for (int i=0; i<3; i++) qvqv += qv(iK,i)*qv(jK,i);
               if (qvs(iK)*qvs(jK)==0) argums(iK,jK)=1.;
               else argums(iK,jK) = qvqv/(qvs(iK)*qvs(jK));
            }
         }
    """
    qv = np.zeros((len(Km),3), dtype=float)
    qvs = np.zeros(len(Km), dtype=float)
    argums = np.zeros((len(Km), len(Km)), dtype=float)
    weave.inline(codeArg, ['Km', 'k', 'qv', 'qvs', 'argums'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    
    # This part of code is too slow in Python even when using special.Legendre
    # It computes the Legendre polynomial for all values of argums(K,K') precomputes above
    # and for all l up to lmax=len(Enu)
    # The few lowest order Legendre polynomials are precomputes and recursion is used only for high order l.
    #   Computes:  Leg(K,K',l) = P_l(argums(K,K'))
    codeLegendre="""
      for (int iK=0; iK<argums.shape()[0]; iK++){
         for (int jK=0; jK<argums.shape()[1]; jK++){
            double x=argums(iK,jK);
            double x2 = x*x;
            Leg(iK,jK,0) = 1;
            if (lmax>=1)  Leg(iK,jK,1) = x;
            if (lmax>=2)  Leg(iK,jK,2) = 1.5*x2-0.5;
            if (lmax>=3)  Leg(iK,jK,3) = x*(2.5*x2-1.5);
            if (lmax>=4)  Leg(iK,jK,4) = 0.375*(1-10*x2*(1-1.1666666666666667*x2));
            if (lmax>=5)  Leg(iK,jK,5) = 1.875*x*(1-4.66666666666666667*x2*(1-0.9*x2));

            for (int l=6; l<=lmax; l++){
                double p0 = 0.375*(1-10*x2*(1-1.1666666666666667*x2));
                double p1 = 1.875*x*(1-4.66666666666666667*x2*(1-0.9*x2)); 
                double p2=0;
                for (int i=6; i<=l; i++){
                  p2 = ((2*i-1)*x*p1-(i-1)*p0)/i;
                  p0=p1;
                  p1=p2;
                }
                Leg(iK,jK,l) = p2;
            }
         }
      }
    """
    lmax = len(Enu)-1
    Leg = np.zeros((len(Km),len(Km),len(Enu)), dtype=float)
    weave.inline(codeLegendre, ['argums', 'lmax', 'Leg'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')

    
    # This part of code is too slow in Python
    # Implements Eq. 46 and Eq. 47 on page 26
    # It computes the Hamiltonian and Overlap in Muffin-Thin and interstitials
    # All necessary arrays were precomputes above.
    codeHO="""
        for (int iK=0; iK<Ham.shape()[0]; iK++){
            for (int jK=0; jK<Ham.shape()[1]; jK++){
                double olapMT=0;
                double hamMT=0;
                for (int l=0; l<Enu.size(); l++){
                    double tC2l = C1(iK,l)*C1(jK,l)*Leg(iK,jK,l);
                    double toop = 1 + omegal(iK,l)*omegal(jK,l)*PP(l);
                    olapMT += tC2l*toop;
                    hamMT  += tC2l*(0.5*(omegal(iK,l)+omegal(jK,l))+toop*Enu(l));
                    C2l(l,iK,jK) = tC2l;
                    C2_1(l,iK,jK) = tC2l * (omegal(iK,l)+omegal(jK,l));
                    C2_2(l,iK,jK) = tC2l * (omegal(iK,l)*omegal(jK,l));
                }
                Olap(iK,jK) = olapMT + Olap_I(iK,jK);
                Ham(iK,jK) = ( 0.25*(qvs(iK)*qvs(iK) + qvs(jK)*qvs(jK)) + VKSi )*Olap_I(iK,jK) + hamMT;
            }
        }
    """
    Olap = np.zeros((len(Km), len(Km)), dtype=float)
    Ham  = np.zeros((len(Km), len(Km)), dtype=float)
    C2l = np.zeros((len(Enu), len(Km), len(Km)), dtype=float)
    C2_1 = np.zeros((len(Enu), len(Km), len(Km)), dtype=float)
    C2_2 = np.zeros((len(Enu), len(Km), len(Km)), dtype=float)
    Enu = np.array(Enu)
    weave.inline(codeHO, ['Olap', 'Ham', 'C2l', 'C2_1', 'C2_2', 'Enu', 'C1', 'Leg', 'omegal', 'qvs', 'PP', 'Olap_I', 'VKSi'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    # Diagonalizes the LAPW Hamiltonian
    #Ek, Ar = symeig.symeig(Ham, Olap, type=1) # symmetric generalized eigenvalue problem
    Ek, Ar = scipy.linalg.eigh(Ham, Olap) # symmetric generalized eigenvalue problem

    
    #print matrix(Ar).T * matrix(Ham) * matrix(Ar)
    
    # Calculation of weights for valence density
    # Implements the weight functions Eqs.58, 59, 60, 61 on page 30.
    Ar = np.matrix(Ar)
    w0 = np.zeros((len(Enu),len(Ar)),dtype=float)
    w1 = np.zeros((len(Enu),len(Ar)),dtype=float)
    w2 = np.zeros((len(Enu),len(Ar)),dtype=float)
    for l in range(len(Enu)):
        tw0 = Ar.T * np.matrix(C2l[l])  * Ar
        tw1 = Ar.T * np.matrix(C2_1[l]) * Ar
        tw2 = Ar.T * np.matrix(C2_2[l]) * Ar
        w0[l,:] = np.array([tw0[p,p] for p in range(len(Ar))])
        w1[l,:] = np.array([tw1[p,p] for p in range(len(Ar))])
        w2[l,:] = np.array([tw2[p,p] for p in range(len(Ar))])
    twi = Ar.T * np.matrix(Olap_I) * Ar
    wi = np.array([twi[p,p] for p in range(len(Ar))])

    return (Ek, Ar, w0, w1, w2, wi)
