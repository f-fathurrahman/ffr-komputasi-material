include("../1d/build_D2_matrix_3pt.jl")
include("../1d/build_D2_matrix_5pt.jl")
include("../1d/build_D2_matrix_7pt.jl")
include("../1d/build_D2_matrix_9pt.jl")

const ⊗ = kron

function build_nabla2_matrix( fdgrid::FD3dGrid; func_1d=build_D2_matrix_3pt )
    
    Nx = fdgrid.Nx
    hx = fdgrid.hx

    Ny = fdgrid.Ny
    hy = fdgrid.hy

    Nz = fdgrid.Nz
    hz = fdgrid.hz

    D2x = func_1d(Nx, hx)
    D2y = func_1d(Ny, hy)
    D2z = func_1d(Nz, hz)

    IIx = speye(Nx)
    IIy = speye(Ny)
    IIz = speye(Nz)

    #∇2 = kron(D2x, speye(Ny)) + kron(speye(Nx), D2y)

    ∇2 = D2x⊗IIy⊗IIz + IIx⊗D2y⊗IIz + IIx⊗IIy⊗D2z 

    return ∇2

end