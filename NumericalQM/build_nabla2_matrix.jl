function build_nabla2_matrix( fdgrid::FD2dGrid; build_D2_matrix::Function=build_D2_matrix_3pt )
    Nx = fdgrid.Nx
    hx = fdgrid.hx
    Ny = fdgrid.Ny
    hy = fdgrid.hy
    D2x = build_D2_matrix(Nx, hx)
    D2y = build_D2_matrix(Ny, hy)
    ∇2 = kron(D2x, speye(Ny)) + kron(speye(Nx), D2y)
    return ∇2
end