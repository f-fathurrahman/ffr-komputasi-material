using Printf
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU
using Random

include("init_FD1d_grid.jl")
include("FD2dGrid.jl")
include("build_D2_matrix.jl")
include("build_nabla2_matrix.jl")
include("supporting_functions.jl")
include("diag_Emin_PCG.jl")
include("ortho_sqrt.jl")
include("ortho_gram_schmidt.jl")

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

function main()

    Random.seed!(1234)

    Nx = 50
    Ny = 50
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )

    ∇2 = build_nabla2_matrix( fdgrid, func_1d=build_D2_matrix_3pt )

    Vpot = pot_harmonic( fdgrid )
    
    Ham = -0.5*∇2 + spdiagm( 0 => Vpot )

    #prec = ilu(-0.5*∇2)
    prec = ilu(Ham)

    # solve for 5 lowest (using `false`) eigenvalues
    Nstates = 5
    Npoints = Nx*Ny
    X = rand(Float64, Npoints, Nstates)
    ortho_sqrt!(X)
    evals = diag_Emin_PCG!( Ham, X, prec, verbose=true )
end

main()



