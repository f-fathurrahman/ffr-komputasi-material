function create_D2_matrix_FD(FDn::Int64, Lx::Float64, Nx::Int64, BCx::Symbol)
    #
    if BCx == :DIRICHLET_BC
        dx = Lx/(Nx-1)
    elseif BCx == :PERIODIC_BC
        dx = Lx/Nx
    else
        error("BCx = $BCx is not supported")
    end
    #
    # Finite difference weights of the second derivative
    w2 = zeros(Float64, FDn + 1)
    for k in 1:FDn
        w2[k+1] = (2*(-1)^(k+1))*(factorial(FDn)^2)/(k*k*factorial(FDn-k)*factorial(FDn+k))
        w2[1] = w2[1] - 2*(1/(k*k))
    end
    #
    n0 = FDn
    # Initial number of non-zeros: including ghost nodes
    tmp_nnzCount = (2*n0 + 1) * Nx
    #
    # Row and column indices and the corresponding non-zero values
    # used to generate sparse matrix DL11 s.t. DL11(I(k),II(k)) = V(k)
    tmp_idx_row = zeros(Int64, tmp_nnzCount)
    tmp_idx_col = zeros(Int64, tmp_nnzCount)
    tmp_vv = zeros(Float64, tmp_nnzCount)
    row_cnt = 1
    cnt = 1
    coef_dxx = 1/dx^2
    # Find non-zero entries that use forward difference
    for ii in 1:Nx
        # diagonal element
        tmp_idx_row[cnt] = row_cnt
        tmp_idx_col[cnt] = ii
        tmp_vv[cnt] = w2[1]*coef_dxx
        cnt += 1 # increment
        # off-diagonal elements
        for q in 1:n0
            # ii + q
            tmp_idx_row[cnt] = row_cnt
            tmp_idx_col[cnt] = ii + q
            tmp_vv[cnt] = w2[q+1]*coef_dxx
            cnt += 1 # increment
            # ii - q
            tmp_idx_row[cnt] = row_cnt
            tmp_idx_col[cnt] = ii - q
            tmp_vv[cnt] = w2[q+1]*coef_dxx
            cnt += 1 # increment
        end
        row_cnt += 1
    end
    #
    if BCx == :DIRICHLET_BC
        # Removing outside domain entries
        idx_inside = (tmp_idx_col .>= 1) .&& (tmp_idx_col .<= Nx)
        idx_row = tmp_idx_row[idx_inside]
        idx_col = tmp_idx_col[idx_inside]
        vv = tmp_vv[idx_inside]
    elseif BCx == :PERIODIC_BC
        is_out_left = tmp_idx_col .< 1
        is_out_right = tmp_idx_col .> Nx # Warning: Assumed influence of only neighboring cells
        idx_row = tmp_idx_row[:] # force copy
        idx_col = mod.(tmp_idx_col .+ (Nx-1), Nx) .+ 1
        vv = tmp_vv[:];
    else
        error("BCx = $BCx is not supported")
    end
    #
    # Build the matrix (dense)
    D2mat = zeros(Float64, Nx, Nx)
    nnzCount = length(vv)
    for ip in 1:nnzCount
        i = idx_row[ip]
        j = idx_col[ip]
        D2mat[i,j] = vv[ip]
    end
    #
    return D2mat
end