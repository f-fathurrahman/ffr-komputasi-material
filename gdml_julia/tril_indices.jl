function tril_indices(N)
    Nelements = (N * (N - 1)) / 2 |> Int64
    idx_rows = zeros(Int64, Nelements)
    idx_cols = zeros(Int64, Nelements)
    ip = 1
    # XXX: Use CartesianIndex?
    for irow in 2:N, icol in 1:N
        if irow > icol
            idx_rows[ip] = irow
            idx_cols[ip] = icol
            ip += 1
        end
    end
    return idx_rows, idx_cols
end

function tril_linear_indices(N)
    Nelements = (N * (N - 1)) / 2 |> Int64
    idx_lin = zeros(Nelements)
    ip = 0
    # XXX: Use CartesianIndex?
    for irow in 2:N, icol in 1:N
        if irow > icol
            idx[ip] = CartesianIndex(irow,icol)
            ip += 1
        end
    end
    return idx
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

2,1
3,1
3,2
4,1
4,2
4,3
=#

