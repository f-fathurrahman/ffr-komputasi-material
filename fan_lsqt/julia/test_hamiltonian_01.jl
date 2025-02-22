using Printf
using SparseArrays

using Infiltrator

include("mod_lsqt_01.jl")

function test_main()
    Nx = 50_000
    Ny = 2
    W = 1.0
    H, V = init_H_and_V(Nx, Ny, W)
    @infiltrate
end
