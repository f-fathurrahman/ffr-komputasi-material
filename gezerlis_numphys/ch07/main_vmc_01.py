from vmc_01 import calc_psi, calc_ekin, calc_epot
import numpy as np

def stats(fs):
    n = fs.size
    fbar = np.sum(fs)/n
    fsq = np.sum(fs**2)/n
    varfbar = (fsq - fbar**2)/(n - 1)
    return fbar, np.sqrt(varfbar)

def vmc(npart, ndim, al, oms, inseed=8735):
    Ncal, nm, th = 10**4, 100, 0.8
    np.random.seed(inseed)
    rolds = np.random.uniform(-1, 1, (npart, ndim))
    psiold = calc_psi(al, oms, rolds)
    iacc, imeas = 0, 0
    eners = np.zeros(Ncal)
    for itot in range(nm*Ncal):
        rnews = rolds+th*np.random.uniform(-1,1,(npart, ndim))
        psinew = calc_psi(al, oms, rnews)
        psiratio = (psinew/psiold)**2
        if psiratio >= np.random.uniform(0,1):
            rolds = np.copy(rnews)
            psiold = psinew
            iacc +=1
        if (itot%nm)==0:
            eners[imeas] = calc_ekin(al,oms,rolds) + calc_epot(oms,rolds)
            imeas += 1
    return iacc/(nm*Ncal), eners

if __name__ == '__main__':
    npart = 4
    ndim = 3
    al = 0.6
    oms = np.arange(1, 1 + ndim)
    accrate, eners = vmc(npart, ndim, al, oms)
    av, err = stats(eners)
    print(accrate, av, err)

