function create_D1_matrix_FD(FDn::Int64, Lx::Float64, Nx::Int64, BCx::Symbol)
    #
    if BCx == :DIRICHLET_BC
        dx = Lx/(Nx-1)
    elseif BCx == :PERIODIC_BC
        dx = Lx/Nx
    else
        error("BCx = $BCx is not supported")
    end
    #
    w1 = zeros(Float64, FDn + 1)
    for k in 1:FDn
        w1[k+1] = ((-1)^(k+1))*(factorial(FDn)^2)/(k*factorial(FDn-k)*factorial(FDn+k))
    end
    #
    # Alias
    n0 = FDn
    # 1st derivative matrix
    nnz_x = 2*n0*Nx
    tmp_idx_row = zeros(Int64, nnz_x)
    tmp_idx_col = zeros(Int64, nnz_x)
    tmp_vv = zeros(Float64, nnz_x)
    row_cnt = 1
    cnt = 1
    coef_dx = 1/dx
    #
    for ii in 1:Nx
        for q in 1:n0
            # ii + q
            tmp_idx_row[cnt] = row_cnt
            tmp_idx_col[cnt] = ii + q
            tmp_vv[cnt] = w1[q+1]*coef_dx
            cnt += 1
            # ii - q
            tmp_idx_row[cnt] = row_cnt
            tmp_idx_col[cnt] = ii - q
            tmp_vv[cnt] = -w1[q+1]*coef_dx # antisymmetric
            cnt += 1
        end
        row_cnt = row_cnt + 1
    end

    if BCx == :DIRICHLET_BC
        # Removing outside domain entries
        idx_inside = (tmp_idx_col .>= 1) .&& (tmp_idx_col .<= Nx)
        idx_row = tmp_idx_row[idx_inside]
        idx_col = tmp_idx_col[idx_inside]
        vv = tmp_vv[idx_inside]
    elseif BCx == :PERIODIC_BC
        is_out_left = tmp_idx_col .< 1
        is_out_right = tmp_idx_col .> Nx
        idx_row = tmp_idx_row[:] # force copy
        idx_col = mod.(tmp_idx_col .+ (Nx-1), Nx) .+ 1
        vv = tmp_vv[:];
    else
        error("BCx = $BCx is not supported")
    end

    # Build the matrix (dense)
    D1mat = zeros(Float64, Nx, Nx)
    nnzCount = length(vv)
    for ip in 1:nnzCount
        i = idx_row[ip]
        j = idx_col[ip]
        D1mat[i,j] = vv[ip]
    end
    #
    return D1mat
    #
end # function