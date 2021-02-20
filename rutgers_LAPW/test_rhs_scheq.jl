using Printf

include("rhs_scheq.jl")
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
end

main()