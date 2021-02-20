using Printf

include("integ_simpson13.jl")

function f(x)
    return 3*cos(x)
end

function main()
    A = 0.0
    B = 1.0
    Npoints = 201
    x = range(A, stop=B, length=Npoints)
    fx = f.(x)

    s = integ_simpson13(fx, x)
    println("s = ", s)
end

main()