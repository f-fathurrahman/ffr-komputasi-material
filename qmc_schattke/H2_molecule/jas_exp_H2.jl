
function jas_exp_H2!(
    global_vars,
    idx_electron,
    vj, vjd, v, vd
)

    # Updates the distances from the active electron to all others
    # and determines Jastrow exponent and derivatives for atoms:
    # u = F/2/(1+r/F)*( delta_{s,-s'} + 1/2*delta_{s,s'} )
    # In addition, the wave function of Kolos and Roothaan
    # (James and Coolidge)
    # (Rev.Mod..Phys.32,219(1960)) is programmed as part of the Jastrow
    # factor. It has to be combined with "product ansatz" of orbital wave.
    
    R_electrons = global_vars.R_electrons
    R_electrons_new = global_vars.R_electrons_new

    a = zeros(Float64, 12)
    # Without series expansion comment 12 lines of constants from KR
    # Constants from Kolos and Roothaan
    a[1]  = +2.192089
    a[2]  = +1.098975
    a[3]  = -0.377500
    a[4]  = -0.139338
    a[5]  = +0.859247
    a[6]  = -0.058316
    a[7]  = +0.078257
    a[8]  = +0.150633
    a[9]  = -0.052156
    a[10] = -0.126629
    a[11] = +0.132561
    a[12] = +0.248411

    jasu = 0.0
    jasup = 0.0
    jasdif = 0.0
    jasdifp = 0.0
    u2d = 0.0
    u2do = 0.0
    u3 = 0.0
    u4 = 0.0

    u1d   = zeros(Float64, 3)
    u1dp  = zeros(Float64, 3)
    u1do  = zeros(Float64, 3)
    u1dpo = zeros(Float64, 3)
    dR = zeros(Float64, 4, Nelectrons, Nelectrons)
    dR_new = zeros(Float64, 4, Nelectrons, Nelectrons)

    EMACH = 1.0e-8 # a quite small number

    for k in 1:Nelectrons
        if k == idx_electron
            continue
        end
        
        as = 0.5     # for equal spins
        
        # This check whether idx_electron and k are of different spin
        cond1 = (idx_electron <= NelectronsPerSpin) && (k > NelectronsPerSpin)
        cond2 = (idx_electron >  NelectronsPerSpin) && (k <= NelectronsPerSpin)
        if (cond1 || cond2)
            as = 1.0
        end
        #
        dR[1:3,idx_electron,k] .= R_electrons[1:3,idx_electron] .- R_electrons[1:3,k]
        dR_new[1:3,idx_electron,k] .= R_electrons_new[1:3,idx_electron] .- R_electrons[1:3,k]
        # Restrict lower limit of the length
        woo = max(2.0/3.0*EMACH, sqrt( sum(dR[1:3,idx_electron,k].^2)) )
        won = max(2.0/3.0*EMACH, sqrt( sum(dR_new[1:3,idx_electron,k].^2)) )
        dR[4,idx_electron,k] = woo
        dR_new[4,idx_electron,k] = won
        #
        jasdif += CJAS*( as/2.0/(1.0 + won/CJAS) - as/2.0/(1.0 + woo/CJAS) ) +
                  β1*( exp(-β2*won^2) - exp(-β2*woo^2) )
    
        jasdifp += sum( ( R_electrons_new[1:3,idx_electron] + R_electrons[1:3,k] - R_atoms[1:3,2] )^2 ) -
                   sum( ( R_electrons[1:3,idx_electron] + R_electrons[1:3,k] - R_atoms[1:3,2] )^2 )
        
        jasu += as*CJAS/2.0/(1.0 + won/CJAS) + β1*exp(-β2*won^2)
    
        jasup += sum( (R_electrons_new[1:3,idx_electron] + R_electrons[1:3,k] - R_atoms[1:3,2])^2 )
    
        u1d[1:3] .-= as*dR_new[1:3,idx_electron,k] / won / (1.0 + won/CJAS)^2/2.0 .-
                     2.0*dR_new[1:3,idx_electron,k]*β1*β2*exp(-β2*won^2)
    
        u1dp[1:3] = R_electrons_new[1:3,idx_electron] + R_electrons[1:3,k] - R_atoms[1:3,2]
        u1do[1:3] .-= as*dR[1:3,idx_electron,k] / woo / (1.0 + woo/CJAS)^2/2.0 .-
                      2.0*dR[1:3,idx_electron,k]*β1*β2*exp(-β2*woo^2)
        u1dpo[1:3] = R_electrons[1:3,idx_electron] + R_electrons[1:3,k] - R_atoms[1:3,2]
    
        u2d -= as/won/(1.0 + won/CJAS)^3 + β1*β2*exp(-β2*won^2)*(4.0*β2*won^2 - 6.0)
  
        u2do -= as/woo/(1.0+woo/CJAS)^3 + β1*β2*exp(-β2*woo^2)*(4.0*β2*woo^2 - 6.0)
        
        u3 += 1.0/won
        u4 += 1.0/won - 1.0/woo
    end
  
    jasu += jasup*γ/2.0
    jasdif += jasdifp*γ/2.0
    u1d += u1dp*γ
    u1do += u1dpo*γ
    u2d += 3.0*γ
    u2do += 3.0*γ
    
    # For an additional Jastrow factor: psi_J=J*sum_{ii=1}^5 a(ii)*x(ii):
    # The basis functions x(ii) are formulated in terms of
    # lambda1,lambda2,mu1,mu2,rho from James and Coolidge, which are abbreviated
    # by the first two letters; the letters a and b denote the two protons and
    # the letters n and o are appended to refer to new and old;
    # the suffixes 1 and 2 are displayed by the
    # index idx_electron and kie of the actual electron which might have been moved
    # and the other second electron which keeps its old position, resp..
    kie = 2
    if idx_electron == 2
        kie = 1
    end
    
    rhn = 2.0*dR_new[4,idx_electron,kie]/DKX
    rho = 2.0*dR[4,idx_electron,kie]/DKX
    
    ran[idx_electron] = max(EMACH, sqrt(sum( (R_electrons_new[1:3,idx_electron] - R_atoms[1:3,1])^2 )))
    rao[idx_electron] = max(EMACH, sqrt(sum( (R_electrons[1:3,idx_electron]   - R_atoms[1:3,1])^2 )))
    rbn[idx_electron] = max(EMACH, sqrt(sum( (R_electrons_new[1:3,idx_electron] - R_atoms[1:3,2])^2 )))
    rbo[idx_electron] = max(EMACH, sqrt(sum( (R_electrons[1:3,idx_electron]   - R_atoms[1:3,2])^2 )))
    
    lan[idx_electron] = ( ran[idx_electron] + rbn[idx_electron] )/DKX
    lao[idx_electron] = ( rao[idx_electron] + rbo[idx_electron] )/DKX
    mun[idx_electron] = ( ran[idx_electron] - rbn[idx_electron] )/DKX
    muo[idx_electron] = ( rao[idx_electron] - rbo[idx_electron] )/DKX
  
    rao[kie] = max(EMACH, sqrt(sum( (R_electrons[1:3,kie] - R_atoms[1:3,1])^2 )))
    rbo[kie] = max(EMACH, sqrt(sum( (R_electrons[1:3,kie] - R_atoms[1:3,2])^2 )))
    lao[kie] = (rao[kie] + rbo[kie])/DKX
    muo[kie] = (rao[kie] - rbo[kie])/DKX
    
    # Accepted step: new coordinates for idx_electron and old for kie
    xn[1] = 2.0
    xn[2] = mun[idx_electron]^2 + muo[kie]^2
    xn[3] = 2.0*mun[idx_electron]*muo[kie]
    xn[4] = lan[idx_electron] + lao[kie]
    xn[5] = 2.0*rhn
    xn[6] = (lan[idx_electron] + lao[kie])*mun[idx_electron]*muo[kie]
    xn[7] = lan[idx_electron]*muo[kie]^2 + lao[kie]*mun[idx_electron]^2
    xn[8] = lao[kie]^2 + lan[idx_electron]^2
    xn[9] = 2.0*rhn^2
    xn[10] = 2.0*lan[idx_electron]*lao[kie]
    xn[11] = 2.0*mun[idx_electron]^2*muo[kie]^2
    xn[12] = (muo[kie]^2 + mun[idx_electron]^2)*rhn
    
    # Not accepted step: old coordinates for both idx_electron and kie
    xo[1] = 2.0
    xo[2] = muo[idx_electron]^2 + muo[kie]^2
    xo[3] = 2.0*muo[idx_electron]*muo[kie]
    xo[4] = lao[idx_electron] + lao[kie]
    xo[5] = 2.0*rho
    xo[6] = (lao[idx_electron] + lao[kie])*muo[idx_electron]*muo[kie]
    xo[7] = lao[idx_electron]*muo[kie]^2 + lao[kie]*muo[idx_electron]^2
    xo[8] = lao[kie]^2 + lao[idx_electron]^2
    xo[9] = 2.0*rho^2
    xo[10] = 2.0*lao[idx_electron]*lao[kie]
    xo[11] = 2.0*muo[idx_electron]^2*muo[kie]^2
    xo[12] = (muo[kie]^2 + muo[idx_electron]^2)*rho
    
    # The 1st derivative (new and old) is denoted by yn(3,12) and yo(3,12)
    @. rran[1:3] = ( R_electrons_new[1:3,idx_electron] - R_atoms[1:3,1])/DKX/ran[idx_electron]
    @. rrbn[1:3] = ( R_electrons_new[1:3,idx_electron] - R_atoms[1:3,2])/DKX/rbn[idx_electron]
    @. rrao[1:3] = ( R_electrons[1:3,idx_electron] - R_atoms[1:3,1])/DKX/rao[idx_electron]
    @. rrbo[1:3] = ( R_electrons[1:3,idx_electron] - R_atoms[1:3,2])/DKX/rbo[idx_electron]
    
    #Accepted step
    yn[1:3,1]  = 0.0
    yn[1:3,2]  = 2.0*mun[idx_electron]*(rran - rrbn)
    yn[1:3,3]  = 2.0*muo[kie]*(rran - rrbn)
    yn[1:3,4]  = rran + rrbn
    yn[1:3,5]  = 4.0*dR_new[1:3,idx_electron,kie]/dR_new[4,idx_electron,kie]/DKX
    yn[1:3,6]  = (rran + rrbn)*mun[idx_electron]*muo[kie] + (lan[idx_electron] + lao[kie])*(rran - rrbn)*muo[kie]
    yn[1:3,7]  = (rran + rrbn)*muo[kie]^2 + 2.0*lao[kie]*mun[idx_electron]*(rran - rrbn)
    yn[1:3,8]  = 2.0*lan[idx_electron]*(rran + rrbn)
    yn[1:3,9]  = 16.0*dR_new[1:3,idx_electron,kie]/DKX^2
    yn[1:3,10] = 2.0*(rran + rrbn)*lao[kie]
    yn[1:3,11] = 4.0*(rran - rrbn)*mun[idx_electron]*muo[kie]^2
  
    yn[1:3,12] = 2.0*mun[idx_electron]*(rran - rrbn)*rhn + 
                 (muo[kie]^2 + mun[idx_electron]^2)*2.0*dR_new[1:3,idx_electron,kie] / dR_new[4,idx_electron,kie]/DKX

    # Not accepted step
    yo[1:3,1] = 0.0
    yo[1:3,2] = 2.0*muo[idx_electron]*(rrao-rrbo)
    yo[1:3,3] = 2.0*muo[kie]*(rrao-rrbo)
    yo[1:3,4] = rrao + rrbo
    yo[1:3,5] = 4.0*dR(1:3,idx_electron,kie)/dR(4,idx_electron,kie)/DKX     
    yo[1:3,6] = (rrao + rrbo)*muo[idx_electron]*muo[kie] + (lao[idx_electron] + lao[kie])*(rrao - rrbo)*muo[kie]
    yo[1:3,7] = (rrao+rrbo)*muo[kie]^2 + 2.0*lao[kie]*muo[idx_electron]*(rrao-rrbo)
    yo[1:3,8] = 2.0*lao[idx_electron]*(rrao + rrbo)
    yo[1:3,9] = 16.0*dR[1:3,idx_electron,kie]/DKX^2
    yo[1:3,10] = 2.0*(rrao + rrbo)*lao[kie]
    yo[1:3,11] = 4.0*(rrao - rrbo)*muo[idx_electron]*muo[kie]^2
  
    yo[1:3,12] = 2.0*muo[idx_electron]*(rrao - rrbo)*rho +
            (muo[kie]^2+muo[idx_electron]^2)*2.0*dR[1:3,idx_electron,kie]/dR[4,idx_electron,kie]/DKX
    
    # The 2nd derivative (new and old) is denoted by zn(5) and zo(5)
    g1mrn = dot( dR_new[1:3,idx_electron,kie], (rran[1:3] - rrbn[1:3]) ) * 4.0/DKX^2/rhn
  
    g1mro = dot( dR[1:3,idx_electron,kie], (rrao[1:3] - rrbo[1:3]) ) * 4.0/DKX^2/rho
  
    g1lmn = dot(rran,rrbn)*DKX^2
    g1lmo = dot(rrao,rrbo)*DKX^2
    
    g2ln = 2.0*(1.0 + g1lmn)/DKX^2
    g2lo = 2.0*(1.0 + g1lmo)/DKX^2
    g2mn = 2.0*(1.0 - g1lmn)/DKX^2
    g2mo = 2.0*(1.0 - g1lmo)/DKX^2
    lapln = 2.0*(1.0/ran[idx_electron] + 1.0/rbn[idx_electron])/DKX
    laplo = 2.0*(1.0/rao[idx_electron] + 1.0/rbo[idx_electron])/DKX
    lapmn = 2.0*(1.0/ran[idx_electron] - 1.0/rbn[idx_electron])/DKX
    lapmo = 2.0*(1.0/rao[idx_electron] - 1.0/rbo[idx_electron])/DKX
  
    # Accepted step
    zn[1] = 0.0
    zn[2] = 2.0*(g2mn + mun[idx_electron]*lapmn)
    zn[3] = 2.0*muo[kie]*lapmn
    zn[4] = lapln
    zn[5] = 16.0/rhn/DKX^2
    zn[6] = muo[kie]*(lapln*mun[idx_electron] + lapmn*(lan[idx_electron] + lao[kie]))
    zn[7] = lapln*muo[kie]^2 + 2.0*lao[kie]*(g2mn + mun[idx_electron]*lapmn)
    zn[8] = 2.0*(g2ln + lan[idx_electron]*lapln)
    zn[9] = 48.0/DKX^2
    zn[10] = 2.0*lao[kie]*lapln
    zn[11] = 4.0*muo[kie]^2*(g2mn + mun[idx_electron]*lapmn)
    zn[12] = 2.0*rhn*(g2mn+mun[idx_electron]*lapmn) + 4.0*mun[idx_electron]*g1mrn + 8.0*(mun[idx_electron]^2 + muo[kie]^2)/rhn/DKX^2

    # Not accepted step
    zo[1] = 0.0
    zo[2] = 2.0*(g2mo+muo[idx_electron]*lapmo)
    zo[3] = 2.0*muo[kie]*lapmo
    zo[4] = laplo
    zo[5] = 16.0/rho/DKX^2
    zo[6] = muo[kie]*(laplo*muo[idx_electron]+lapmo*(lao[idx_electron]+lao[kie]))
    zo[7] = laplo*muo[kie]^2+2.0*lao[kie]*(g2mo+muo[idx_electron]*lapmo)
    zo[8] = 2.0*(g2lo+lao[idx_electron]*laplo)
    zo[9] = 48.0/DKX^2
    zo[10] = 2.0*lao[kie]*laplo
    zo[11] = 4.0*muo[kie]^2*(g2mo+muo[idx_electron]*lapmo)
    zo[12] = 2.0*rho*(g2mo + muo[idx_electron]*lapmo) + 4.0*muo[idx_electron]*g1mro + 8.0*(muo[idx_electron]^2+muo[kie]^2)/rho/DKX^2
    
    jasjcn = 0.0
    jasjco = 0.0
    jc1dn[1:3] = 0.0
    jc1do[1:3] = 0.0
    jc2dn = 0.0
    jc2do = 0.0
  
    for ii in 1:12
        jasjcn += a[ii]*xn[ii]
        jasjco += a[ii]*xo[ii]
        jc1dn[1:3] .+= a[ii]*yn[1:3,ii]
        jc1do[1:3] .+= a[ii]*yo[1:3,ii]
        jc2dn += a[ii]*zn[ii]
        jc2do += a[ii]*zo[ii]
    end
  
    # if((jasjcn .le. 0.0) .or. (jasjco .le. 0.0)) then
    #   write(*,*)'non-positive argument of log'
    #   stop
    # endif

    # jasjc or jasjco might be negative, but only the square is needed
    # for use in probability measure. So, take modulus.
    jasujcn = -log(abs(jasjcn))
    jasujco = -log(abs(jasjco))
    
    #Instead calculate directly the acceptance ratio qjc
    QJC = jasjcn/jasjco 

    # For derivatives wave function is just a factor, no exponentiation
    jc1dn[1:3] = jc1dn[1:3]/jasjcn
    jc1do[1:3] = jc1do[1:3]/jasjco
    jc2dn = jc2dn/jasjcn
    jc2do = jc2do/jasjco
    vjd[idx_electron] = jasdif + jasujcn - jasujco
    vj[idx_electron] = jasu + jasujcn
    GRJAS[1:3,idx_electron] = - u1d[1:3] + jc1dn[1:3]
    LAPJAS[idx_electron] = -u2d + sum(u1d[1:3]^2) + jc2dn - 2.0*dot(jc1dn[1:3],u1d[1:3])
    GRJASOLD[1:3,idx_electron] = - u1do[1:3] + jc1do[1:3]

    LAPJASOLD[idx_electron] = -u2do + sum(u1do[1:3]^2) + jc2do - 2.0*dot( jc1do[1:3], u1do[1:3] )
    v[idx_electron] = u3
    vd[idx_electron] = u4

    return
end
