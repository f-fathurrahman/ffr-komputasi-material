# using column major
function idx_rowcol_to_linear(i,j,N)
    return i + (j-1)*N
end


function tril_indices(N)
    Nelements = (N * (N - 1)) / 2 |> Int64
    idx_rows = zeros(Int64, Nelements)
    idx_cols = zeros(Int64, Nelements)
    idx_lin = zeros(Int64, Nelements)
    ip = 1
    for icol in 1:N, irow in (icol+1):N
        idx_rows[ip] = irow
        idx_cols[ip] = icol
        idx_lin[ip] = idx_rowcol_to_linear(irow, icol, N)
        ip += 1
    end
    return idx_rows, idx_cols, idx_lin
end


function tril_cart_indices(N)
    Nelements = (N * (N - 1)) / 2 |> Int64
    idx = Vector{CartesianIndex}(undef,Nelements)
    ip = 1
    # XXX: Use CartesianIndex?
    for irow in 2:N, icol in 1:N
        if irow > icol
            idx[ip] = CartesianIndex(irow,icol)
            ip += 1
        end
    end
    return idx
end



#=
0 0 0 0    
1 0 0 0   
1 1 0 0
1 1 1 0

 1  5   9  13
 2  6  10  14
 3  7  11  15
 4  8  12  16

2,1 -> 2 = 2 + (1-1)*4
3,1 -> 3 = 3 + (1-1)*4
3,2 -> 7 = 3 + (2-1)*4
4,1 -> 4 = 4 + (1-1)*4
4,2 -> 8 = 4 + (2-1)*4
4,3 -> 12 = 4 + (3-1)*4
=#

