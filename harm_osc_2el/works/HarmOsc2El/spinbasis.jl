struct SpinBasis{T} <: SpatialBasis
    base::T # The basis with no spin
    l::Int64 # number of basis (2 times larger than spatial basis)
    function SpinBasis(base::SpatialBasis)
        @info "Constructing SpinBasis"
        return new{typeof(base)}(base, 2 * base.l)
    end
end

function spatial(basis::SpinBasis, grid)
    n = length(grid)
    l = basis.l
    res = [zeros(n) for i in 1:l]
    
    nospin = zeros(l÷2)
    
    # loop over spatial points
    for i in 1:n
        #
        evaluate!(nospin, grid[i], basis.base) # the basis functions evaluated at x
        for j in 1:l÷2
            
            #if j % 4 == 0 || (j + 1)%4 == 0 # In case I need the same spfs as ODQD from quantum_systems
            #    nospin[j] = -nospin[j]
            #end
            
            res[2j-1][i] = nospin[j]
            res[2j][i] = nospin[j]
        end
    end
    
    return res
end

