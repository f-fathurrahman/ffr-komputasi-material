# Uncompressed R_d
function uncompress_R_d(Natoms, R_d_desc)
    desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64
    tmp_R_d = zeros(Float64, 3, Natoms, desc_dim)
    idx_rows, idx_cols, idx_lin = tril_indices(Natoms)
    for ip in 1:desc_dim
        i = idx_rows[ip]
        j = idx_cols[ip]
        @views tmp_R_d[:,i,ip] .=  R_d_desc[:,ip]
        @views tmp_R_d[:,j,ip] .= -R_d_desc[:,ip]
    end
    R_d_full = reshape(tmp_R_d, 3*Natoms, desc_dim)
    # Use SparseArray ???
    return R_d_full
end

# use a precalculated TrilIndices
function uncompress_R_d(Natoms, indices::TrilIndices, R_d_desc)
    desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64
    tmp_R_d = zeros(Float64, 3, Natoms, desc_dim)
    #
    idx_rows = indices.idx_rows
    idx_cols = indices.idx_cols
    idx_lin = indices.idx_lin
    #
    for ip in 1:desc_dim
        i = idx_rows[ip]
        j = idx_cols[ip]
        @views tmp_R_d[:,i,ip] .=  R_d_desc[:,ip]
        @views tmp_R_d[:,j,ip] .= -R_d_desc[:,ip]
    end
    R_d_full = reshape(tmp_R_d, 3*Natoms, desc_dim)
    # Use SparseArray ???
    return R_d_full
end
# TODO: use preallocated R_d_full, no need to reshape at the end?
