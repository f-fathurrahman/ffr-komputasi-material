using Printf

include("rhs_scheq.jl")
include("integ_numerov.jl")

function main()
    Nrgrid = 3001
    RmaxAtom = 10.0
    rgrid = collect(range(1e-10, stop=RmaxAtom, length=Nrgrid))
    
    Veff = -1.0./rgrid
    rhs = zeros(Float64,Nrgrid)

    E = -1.0
    l = 1
    rhs_scheq!(E, l, rgrid, Veff, rhs)

    for ir in 1:4
        @printf("%5d %18.10f %18.10e %18.10e\n", ir, rgrid[ir], Veff[ir], rhs[ir])
    end
    ir = Nrgrid-1
    @printf("%5d %18.10f %18.10e %18.10e\n", ir, rgrid[ir], Veff[ir], rhs[ir])
    ir = Nrgrid
    @printf("%5d %18.10f %18.10e %18.10e\n", ir, rgrid[ir], Veff[ir], rhs[ir])

    u = zeros(Float64, Nrgrid)
    h = rgrid[2] - rgrid[1]
    f1 = rgrid[1]*exp(-rgrid[1])
    f2 = rgrid[2]*exp(-rgrid[2])
    #
    integ_numerov!(rhs, h, f1, f2, u)
    #
    for ir in 1:4
        @printf("%5d %18.10e\n", ir, u[ir])
    end
    ir = Nrgrid-1
    @printf("%5d %18.10e\n", ir, u[ir])
    ir = Nrgrid
    @printf("%5d %18.10e\n", ir, u[ir])
end


main()