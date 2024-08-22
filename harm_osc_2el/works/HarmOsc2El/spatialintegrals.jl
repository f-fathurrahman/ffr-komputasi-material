function onebody(ho::HOBasis, grid)
    (; l, ω) = ho
    return LinearAlgebra.Diagonal([(n + 0.5) * ω for n in 0:l-1])
end

function twobody(basis::SpatialBasis, grid, V::Interaction)
    spfs = spatial(basis, grid)
    inner = inner_ints(spfs, grid, V)
    u = outer_int(spfs, grid, inner)
    return u
end

function outer_int(spfs, grid, inner_ints)
    l = length(spfs)
    outer_int = zeros(typeof(spfs[1][1]), l, l, l, l)
    
    fs = [similar(spfs[1]) for i in 1:Threads.nthreads()]
    is = [similar(spfs[1]) for i in 1:Threads.nthreads()]

    @inbounds Threads.@threads for κ in 1:l
        f_vals = fs[Threads.threadid()]
        inner = is[Threads.threadid()]
        
        for λ in 1:l
            @views inner .= inner_ints[κ, λ, :]
            for μ in 1:l
                for ν in 1:l
                    f_vals .= conj.(spfs[μ]) .* inner .* spfs[ν]
                    outer_int[μ, κ, ν, λ] = trapz(f_vals, grid)
                end
            end
        end
    end
    return outer_int
end

function inner_ints(spfs, grid, V::Interaction)
    l = length(spfs)
    inner_int = zeros(typeof(spfs[1][1]), l, l, length(spfs[1]))
    interactions = [similar(grid) for i in 1:Threads.nthreads()]
    fs = [similar(spfs[1]) for i in 1:Threads.nthreads()]
    @inbounds Threads.@threads for xi in eachindex(grid)
        x1 = grid[xi]
        f_vals = fs[Threads.threadid()] # Pre-allocated vector for this thread
        interaction = interactions[Threads.threadid()] # ^^^
        
        interaction = interaction_over_grid!(interaction, x1, grid, V)
        for κ in 1:l
            for λ in κ:l
                f_vals .= conj.(spfs[κ]) .* interaction .* spfs[λ]
                res = trapz(f_vals, grid)
                inner_int[κ, λ, xi] = res
                inner_int[λ, κ, xi] = res'
            end
        end
    end
    return inner_int
end

function inner_ints(spfs, grid, V::NonInteracting)
    l = length(spfs)
    inner_int = zeros(typeof(spfs[1][1]), l, l, length(spfs[1]))
    return inner_int
end

function trapz(f_vals, grid)
    val = sum(f_vals)
    val = val - 0.5 * (f_vals[1] + f_vals[end])
    return val * (grid[2] - grid[1])
end