include("tril_indices.jl")

N = 5
idx_row, idx_col, idx_lin = tril_indices_v2(N)
A = zeros(Float64, N, N)
A[idx_lin] .= 1

B = zeros(Float64, N, N)
for (i,j) in zip(idx_row,idx_col)
    B[i,j] = 7.0
end

