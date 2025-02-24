function my_CSC_mul!(Ax, A, x)
    N = A.n
    @assert N == A.m

    @assert size(x, 1) == N

    # Initialize output vector to zero
    for i in 1:N
        Ax[i] = zero(eltype(A))
    end
    
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval
    # Loop over columns
    for j in 1:N
        for i in colptr[j]:(colptr[j+1]-1)
            Ax[rowval[i]] += nzval[i] * x[j]
        end
    end
    return
end