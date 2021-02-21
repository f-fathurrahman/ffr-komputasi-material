function integ_numerov!(F, dx, f1, f2, solution)
    #
    N = size(F,1)
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
    for  i in 3:N
        w2 = 2*w1 - w0 + h2*phi*Fx
        w0 = w1
        w1 = w2
        Fx = F[i]
        phi = w2/(1.0 - h12*Fx)
        solution[i] = phi
    end
    return
end

# Not yet working
function integ_numerov_inward!(F, dx, fN, fNm1, solution)
    #
    N = size(F,1)
    solution[N] = fN
    solution[N-1] = fNm1
    #
    h2 = -dx*dx
    h12 = h2/12
    #
    w0 = ( 1.0 - h12*F[N] ) * solution[N]
    Fx = F[N-1]
    w1 = ( 1.0 - h12*Fx ) * solution[N-1]
    phi = solution[2]
    for i in (N-2):-1:1
        w2 = 2*w1 - w0 + h2*phi*Fx
        w0 = w1
        w1 = w2
        Fx = F[i]
        phi = w2/(1.0 - h12*Fx)
        solution[i] = phi
    end
    return
end


