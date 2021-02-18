from math import pi

def extrapolate(x, x0, x1, f0, f1):
    "linear extrapolation"
    return f0 + (f1-f0)*(x-x0)/(x1-x0)

def calc_rs(rh):
    "rs from density -> an electron radius that corresponds to density"
    return (3/(4*pi*rh))**(1/3.)