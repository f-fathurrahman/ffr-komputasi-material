using Printf

using PyPlot
const plt = PyPlot

include("rhs_scheq.jl")
include("integ_numerov.jl")

function main()
    Nrgrid = 3001
    RmaxAtom = 10.0
    rgrid = collect(range(1e-5, stop=RmaxAtom, length=Nrgrid))
    
    Veff = -1.0./rgrid
    rhs = zeros(Float64,Nrgrid)

    E = -0.5
    l = 0
    rhs_scheq!(E, l, rgrid, Veff, rhs)

    #for ir in 1:4
    #    @printf("%5d %18.10f %18.10e %18.10e\n", ir, rgrid[ir], Veff[ir], rhs[ir])
    #end
    #ir = Nrgrid-1
    #@printf("%5d %18.10f %18.10e %18.10e\n", ir, rgrid[ir], Veff[ir], rhs[ir])
    #ir = Nrgrid
    #@printf("%5d %18.10f %18.10e %18.10e\n", ir, rgrid[ir], Veff[ir], rhs[ir])

    u = zeros(Float64, Nrgrid)
    h = rgrid[2] - rgrid[1]
    fN = rgrid[Nrgrid]*exp(-rgrid[Nrgrid])
    fNm1 = rgrid[Nrgrid-1]*exp(-rgrid[Nrgrid-1])
    #
    integ_numerov_inward!(rhs, h, fN, fNm1, u)
    #
    for ir in 1:4
        @printf("%5d %18.10e\n", ir, u[ir])
    end
    ir = Nrgrid-2
    @printf("%5d %18.10e\n", ir, u[ir])
    ir = Nrgrid-1
    @printf("%5d %18.10e\n", ir, u[ir])
    ir = Nrgrid
    @printf("%5d %18.10e\n", ir, u[ir])

    plt.clf()
    plt.plot(rgrid, u)
    plt.grid(true)
    plt.savefig("IMG_u_inward.pdf")

    plt.clf()
    plt.plot(rgrid, rhs)
    plt.xlim(0.0, 0.1)
    plt.grid(true)
    plt.savefig("IMG_rhs.pdf")

    println(rgrid[1:5])
    println(rhs[1:5])

end


main()