using LinearAlgebra
using SparseArrays
using BenchmarkTools

include("init_FD1d_grid.jl")
include("FD2dGrid.jl")
include("build_D2_matrix.jl")
include("build_nabla2_matrix.jl")
include("supporting_functions.jl")

function pot_harmonic( fdgrid::FD2dGrid; ω=1.0 )
    Npoints = fdgrid.Npoints
    Vpot = zeros(Npoints)
    for i in 1:Npoints
        x = fdgrid.r[1,i]
        y = fdgrid.r[2,i]
        Vpot[i] = 0.5 * ω^2 *( x^2 + y^2 )
    end
    return Vpot
end

function test_multiplication()
    Nx = 80
    Ny = 80
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    println("Nx = ", Nx)
    println("Ny = ", Ny)

    ∇2 = build_nabla2_matrix( fdgrid )

    psi1 = rand( fdgrid.Nx, fdgrid.Ny )
    psi2 = rand( fdgrid.Npoints )

    println("Using views:")
    @btime begin
       @views res1 = $∇2*$psi1[:]
    end

    println("Not using views:")
    @btime begin
       res1 = $∇2*$psi1[:]
    end

    println("Using one-dimensional array")
    @btime begin
       res2 = $∇2*$psi2
    end



    println("Pass here")
end
test_multiplication()



