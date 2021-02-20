from math import exp

def calc_IS_charge(mu, Ek, wi, wkp, beta=50.):
    " Interstitial charge Eq.64 page 31"
    def ferm(x):
        if x>100: return 0.0
        if x<-100: return 1.0
        return 1/(exp(x)+1.)

    sIntRho=0  #  total interstitial charge
    for ik in range(len(Ek)):
        dsum = 0.0
        for p in range(len(Ek[ik])):
            dsum += ferm( (Ek[ik,p]-mu)*beta )*wi[ik][p]
        sIntRho += dsum*wkp[ik]
    sIntRho *= 2 # due to spin
    return sIntRho
