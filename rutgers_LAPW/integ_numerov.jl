function integ_numerov!(F, dx, f1, f2, solution)
    #
    Nmax = size(F,1)
    solution[1] = f1
    solution[2] = f2
    #
    h2 = dx*dx
    h12 = h2/12
    #
    w0 = ( 1.0 - h12*F[1] ) * solution[1]
    Fx = F[2]
    w1 = ( 1.0 - h12*Fx ) * solution[2]
    phi = solution[2]
    for  i in 3:Nmax
        w2 = 2*w1 - w0 + h2*phi*Fx
        w0 = w1
        w1 = w2
        Fx = F[i]
        phi = w2/(1.0 - h12*Fx)
        solution[i] = phi
    end
    return
end
