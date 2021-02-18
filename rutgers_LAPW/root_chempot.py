import weave

def root_chempot(mu, Ek, wkp, Zval, beta=50.):
    " Computes valence density to find root for the chemical potential "
    code_mu="""
    double Zt=0;
    for (int ik=0; ik<Ek.shape()[0]; ik++){
        for (int p=0; p<Ek.shape()[1]; p++){
            double x = beta*(Ek(ik,p)-mu);
            double ferm = abs(x) < 100 ? 1/(exp(x)+1) : (x<0 ? 1 : 0);
            Zt += wkp(ik) * ferm;
        }
    }
    return_val = Zt;
    """
    Zt = weave.inline(code_mu, ['mu', 'Ek', 'beta', 'wkp'],type_converters=weave.converters.blitz, compiler = 'gcc')
    return 2*Zt-Zval
