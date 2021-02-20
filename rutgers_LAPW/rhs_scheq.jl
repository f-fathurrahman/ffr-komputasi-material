function rhs_scheq!(E, l::Int64, rgrid, Veff, rhs)
    N = size(rgrid,1)
    for i in 1:N
        r = rgrid[i]
        rhs[i] = 2*( -E + 0.5*l*(l+1)/(r*r) + Veff[i] )
    end
    return
end