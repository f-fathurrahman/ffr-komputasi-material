function assemble_Kmatrix!(K, jrow::Int64, Natoms, R_desc_v, R_d_desc_v)

    σ = 20     # kernel parameter, should the the same as task["sig"] ???
    sqrt5 = sqrt(5)
    mat52_base_div = 3*σ^4
    sig_pow2 = σ^2

    Ntrain = length(R_desc_v)

    Rj_d_desc_full = uncompress_R_d(Natoms, R_d_desc_v[jrow])
    
    dim_i = 3*Natoms
    Kij = zeros(Float64, dim_i, dim_i)

    idx_col_start = (jrow-1)*dim_i + 1
    idx_col_stop = jrow*dim_i
    idx_cols = idx_col_start:idx_col_stop

    for irow in jrow:Ntrain

        idx_row_start = (irow-1)*dim_i + 1
        idx_row_stop = irow*dim_i
        idx_rows = idx_row_start:idx_row_stop

        diff_ab = R_desc_v[irow] - R_desc_v[jrow]
        norm_ab = sqrt5 * norm(diff_ab)
        mat52_base = exp(-norm_ab/σ) / mat52_base_div * 5
    
        res1 = diff_ab * mat52_base * 5
        res2 = Rj_d_desc_full * diff_ab # matmul
        diff_ab_outer = res1 * res2'  # matmul
    
        res4 = (sig_pow2 + σ * norm_ab) * mat52_base  # scalar
        diff_ab_outer .-= res4 * Rj_d_desc_full'   # matmul
        Ri_desc_full = uncompress_R_d(Natoms, R_d_desc_v[irow])

        @views Kij[:,:] .= Ri_desc_full * diff_ab_outer
        @views K[idx_rows,idx_cols] .= Kij[:,:]
        @views K[idx_cols,idx_rows] .= Kij[:,:]'
    end

    return

end