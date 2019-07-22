function init_FD1d_grid( X::Tuple{Float64,Float64}, N::Int64 )
    x_min = X[1]
    x_max = X[2]
    L = x_max - x_min
    h = L/(N-1)
    x = zeros(Float64,N)
    for i = 1:N
        x[i] = x_min + (i-1)*h
    end
    return x, h
end