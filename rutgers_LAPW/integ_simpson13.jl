#function integ_simpson13( f, a, b )
#    h = (b - a)/3
#    x0 = a
#    x1 = a + h
#    x2 = b
#    I = h/3 * ( f(x0) + 4*f(x1) + f(x2) )
#    return I
#end
#h/3*(f[1] + 4*f[2] + f[3])
#h/3*(f[3] + 4*f[4] + f[5])
# -> h/3( f[1] + 4*f[2] + 2*f[3] + 4*f[4] + f[5] )

# Assumption: regularly spaced points
function integ_simpson13(f, x)
    N = size(f,1)
    @assert N%2 != 0
    h = x[2] - x[1]
    s = f[1]
    for i in 2:2:N-1
        s = s + 4*f[i]
    end
    for i in 3:2:N-2
        s = s + 2*f[i]
    end
    s = s + f[N]
    return h*s/3
end
