using Printf
using SparseArrays

using Infiltrator

function gen_idx_xy2ip(Nx, Ny)
    idx_xy2ip = zeros(Int64, Nx, Ny)
    ip = 0
    for j in 1:Ny, i in 1:Nx
        ip += 1
        idx_xy2ip[i,j] = ip
    end
    return idx_xy2ip
end


function test_main()

    Nx = 10
    Ny = 10
    W = 1.0


    idx_xy2ip = gen_idx_xy2ip(Nx, Ny)

    Npoints = Nx * Ny # total number of sites

    row_H = zeros(Int64, Npoints*4 - Nx*2) # the row indices for H
    col_H = zeros(Int64, Npoints*4 - Nx*2) # the column indices for H
    Hij = -ones(ComplexF64, Npoints*4 - Nx*2) # nonzero Hamiltonian elements

    row_V = zeros(Int64, Npoints*2) # row indices for V
    col_V = zeros(Int64, Npoints*2) # column indices for V
    Vij = zeros(ComplexF64, Npoints*2) # nonzero velocity matrix elements

    row_U = 1:Npoints # row and column indices for U
    col_U = 1:Npoints
    Uij = 1.234*ones(Npoints) #np.random.uniform(-W * 0.5, W * 0.5, N) # on-site potentials

    count_H = 0
    count_V = 0

    for j in 1:Ny, i in 1:Nx
        #
        # (0) # get the index of the center site
        idx_center = idx_xy2ip[i,j]
        #
        # (1) consider the left neighbor (periodic boundary)
        ii = i - 1
        if ii <= 0
            ii += Nx # wrap
        end
        idx_left = idx_xy2ip[ii,j]
        #
        count_H += 1    
        row_H[count_H] = idx_center
        col_H[count_H] = idx_left
        #
        count_V += 1
        row_V[count_V] = idx_center
        col_V[count_V] = idx_left
        Vij[count_V] = im
        #
        # (2) consider the right neighbor (periodic boundary)
        ii = i + 1
        if ii >= Nx + 1
            ii -= Nx # wrap
        end
        idx_right = idx_xy2ip[ii,j]
        #
        count_H += 1
        row_H[count_H] = idx_center
        col_H[count_H] = idx_right
        #
        count_V += 1
        row_V[count_V] = idx_center
        col_V[count_V] = idx_right
        Vij[count_V] = -im
        #
        # (3) consider the upper neighbor (open boundary)
        @printf("\n -- i=%3d j=%3d idx_center=%d\n", i, j, idx_center)
        @printf("Periodic BC: idx_left=%3d idx_right=%3d\n" , idx_left, idx_right)
        if j < Ny
            idx_up = idx_xy2ip[i,j+1]
            count_H += 1
            row_H[count_H] = idx_center
            col_H[count_H] = idx_up
            @printf("Open BC: idx_up=%3d\n", idx_up)
        end
        #
        # (4) consider the down neighbor (open boundary)
        if j > 1
            idx_down = idx_xy2ip[i,j-1]
            count_H += 1
            row_H[count_H] = idx_center
            col_H[count_H] = idx_down
            @printf("Open BC: index_down=%3d\n", idx_down)
        end
    end

    H = sparse( row_H, col_H, Hij )
    U = sparse( row_U, col_U, Uij )
    H = H + U
    V = sparse( row_V, col_V, Vij )

    @infiltrate

end # function

test_main()
