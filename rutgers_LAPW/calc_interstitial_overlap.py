import numpy as np
from math import pi, sqrt
from scipy import special

def calc_interstitial_overlap(Km, RMuffinTin, Vol):
    """ Overlap in the interstitials can be calculated outside the k-loop
    
        Please see Eq.46 on page 26 for the quantity O_{K'K}^I
    """
    Olap_I = np.zeros((len(Km),len(Km)), dtype=float)
    for i in range(len(Km)):
        Olap_I[i,i] = 1 - 4*pi*RMuffinTin**3/(3.*Vol)
        for j in range(i+1, len(Km)):
            KKl = sqrt(np.dot(Km[i]-Km[j],Km[i]-Km[j]))
            fbessel = special.sph_jn(1,KKl*RMuffinTin)[0][1]
            Olap_I[i,j] = -4*pi*RMuffinTin**2*fbessel/(KKl*Vol)
            Olap_I[j,i] = Olap_I[i,j]

    return Olap_I
