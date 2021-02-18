import numpy as np
from scipy import interpolate

from FccLattice import *
from excor import ExchangeCorrelation
from calc_atom_charge import *
from calc_IS_overlap import *
from solve_scheq import *
from my_utilities import calc_rs
from solve_eigensystem import *
from root_chempot import *
from calc_MT_density import *
from calc_IS_charge import *

DEFAULT_COLOR = '\033[0m'
RED = '\033[31;1m'
GREEN = '\033[32;1m'
BLUE = '\033[34;1m'
YELLOW = '\033[33;1m'

#
# Main start here ...
#

###################################
# Start input parameters
Z = 29                     # Number of electrons in the atom
LatConst = 6.8219117     # Lattic constant
core = [3,2,0]  # 1s,2s,3s, 1p,2p, no-d

nkp = 6                  # Number of k-points in 1BZ: (nkp x nkp x nkp)

#### Core states ##############
#### Linearization energies ###
# Most of high energy partial waves should be centered around mu
# In general, it is a good idea to determine them self-consistently
Enu = [0.11682, 0.18794, 0.211145, 0.3, 0.3, 0.3]
N = 1001                 # Number of points in radial mesh
beta=50.                 # To compute chemical potential we take finite inverse temperature
mu_mm = [0.0, 1.0]       # Chemical potential is searched between mu_mm[0] and mu_mm[1]
CutOffK=3.5              # Largest lengt of reciprocal vectors K (only shorter vec. are taken into account)
DRho = 1e-3              # Convergence criteria for electronic density Rho
Nitt = 100               # Maximum number of itterations
mixRho = 0.3             # Linear mixing parameter for charge
Nkplot = 200             # Number of k-points for plotting bands
plotMM = [-1.,0.1]       # Bands will be plotted in this energy range
# End inout parameters 
########################
    
# Core number of electrons
Zcor = np.sum([2*(2*l+1)*nl for l,nl in enumerate(core)])
# Valence number of electrons
Zval = Z-Zcor
    
print("Z core=", Zcor, " and Zval=", Zval)

# Creates atomic charge to have a good starting point
(Atom_R0, Atom_rho) = calc_atom_charge(Z, core, 0.3)
AtomRhoSpline = interpolate.splrep(Atom_R0, Atom_rho, s=0)

# Exchange correlations class; WVN seems to be the best
# (look http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html)    
XC = ExchangeCorrelation(3)
    
# Generates and stores momentum points
fcc = FccLattice(LatConst)                  # Information about lattice
RMuffinTin = fcc.RMuffinTin()               # Muffin-Tin radius choosen such that spheres touch
VMT = 4*pi*RMuffinTin**3/3.                 # Volume of MT
Vinter = fcc.Volume-VMT                     # Volume of the interstitial region
print("Muffin-Tin radius =", RMuffinTin)
print("Volume of the MT sphere    =", VMT)
print("Volume of the unit cell    =", fcc.Volume)
print("Volume of the interstitial =", Vinter)
fcc.GenerateReciprocalVectors(4, CutOffK)
# Reciprocal bravais lattice is builded, K points taken into account only for |K|<CutOff
fcc.ChoosePointsInFBZ(nkp,0)
# Chooses the path in the 1BZ or the k-points in the irreducible 1BZ

# Radial mesh --  only linear mesh can be used in connection to Numerov algorithm.
R0 = np.linspace(0, RMuffinTin, N)
R0[0]=1e-10
R = R0[::-1]
    
# Interstital overlap does not change through iterations
Olap_I = calc_IS_overlap(fcc.Km, RMuffinTin, fcc.Volume)

# We interpolate atomic charge on the new mesh within Muffin-Tin sphere
TotRho = interpolate.splev(R0, AtomRhoSpline)

for itt in range(Nitt):  # self-consistent loop
    
    print('%d) Preparing potential' % itt)
    UHartree = solve_poisson(Z, R0, TotRho)
    # Adding exchange-correlation part
    Vxc = [XC.Vx(rsi)+XC.Vc(rsi) for rsi in calc_rs(TotRho)]
    #
    nVeff = (UHartree - Z)/R0 + Vxc
    zeroMT = nVeff[-1]  # New MT zero
    nVeff -= zeroMT
    print('   Muffin-Tin zero is ', zeroMT)
    #
    Veff = nVeff
    #
    (logDer, Psi_l, Psip_l) = solve_scheq(Z, Enu, R0, RMuffinTin, Veff)
    #    
    (coreRho, coreE, coreZ, core_states) = find_core_states(core, R0[::-1], Veff[::-1], Z)
    print('   coreZ=', coreZ, 'coreE=', coreE)

    #print("Pass here 103 in main_lapw_py2")
    #exit()
        
    # This is the main loop over all k-points
    Ek=[]; w0=[]; w1=[]; w2=[]; wi=[]
    for ik,k in enumerate(fcc.kp):
        (tEk, tAr, tw0, tw1, tw2, twi) = solve_eigensystem(k, fcc.Km, Olap_I, Enu, logDer, RMuffinTin, fcc.Volume)
        Ek.append(tEk); w0.append(tw0); w1.append(tw1); w2.append(tw2); wi.append(twi);

    Ek = np.array(Ek)

    #print("Pass here 116 in main_lapw_py2")
    #exit()

    # New chemical potential
    mu = scipy.optimize.brentq(root_chempot, mu_mm[0], mu_mm[1], args=(Ek, fcc.wkp, Zval, beta))
    print(GREEN, 'New chemical potential is', mu, DEFAULT_COLOR)
    #
    MTRho = calc_MT_density(mu, Ek, fcc.wkp, w0, w1, w2, Psi_l, Psip_l, beta)
    nTotRho = MTRho + coreRho
    #
    sMTRho = scipy.integrate.simps(MTRho*R0**2*(4*pi), R0)
    sIntRho = calc_IS_charge(mu, Ek, wi, fcc.wkp, beta)
    sCoreRho = scipy.integrate.simps(coreRho*R0**2*(4*pi), R0)
    #    
    print('   Zval=', Zval, '~', sMTRho+sIntRho)
    print('   Weght in the MT sphere =', sMTRho, 'and in the interstitials =', sIntRho, 'and in core =', sCoreRho)
        
    renorm = Z/(sMTRho+sIntRho+sCoreRho)
    print('   Total charge found=', sMTRho+sIntRho+sCoreRho, 'should be', Z, '-> renormalizing by', renorm)
    nTotRho *= renorm

    DiffRho = scipy.integrate.simps(abs(nTotRho-TotRho), R0)
    print(BLUE, 'Electron density difference=', DiffRho, DEFAULT_COLOR)
    if (DiffRho<DRho): break
        
    TotRho = mixRho * nTotRho + (1-mixRho)*TotRho
    


"""
# Plotting bands
fcc.ChoosePointsInFBZ(Nkplot, type=1) 
    
Ek=[]
for ik,k in enumerate(fcc.kp):
    (tEk, tAr, tw0, tw1, tw2, twi) = ComputeEigensystem(k, fcc.Km, Olap_I, Enu, logDer, RMuffinTin, fcc.Volume)
    Ek.append(tEk)
#
Ek = array(Ek)


for i in range(np.shape(Ek)[1]):
    if max(Ek[:,i])-mu > plotMM[0] and min(Ek[:,i])-mu < plotMM[1]:
        plot(Ek[:,i]-mu, 'k-', lw=2)

plot([0,len(Ek)],[0,0], 'k:')  # chemical potential line
ax=axis()

xlabs = [p[1] for p in fcc.Points]
labs  = [p[0] for p in fcc.Points]
xticks(xlabs, labs)
    
for ix,x in enumerate(xlabs):
    plot([x,x], [ax[2],ax[3]], 'k:')
        
axis([xlabs[0], xlabs[-1], ax[2], ax[3]])
show()
"""