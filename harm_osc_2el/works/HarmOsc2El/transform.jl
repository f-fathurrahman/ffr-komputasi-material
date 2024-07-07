function System(system::SpatialSystem, C)
    (; n, l, h, u, spfs, grid, basis, transform, V) = system
    spfs = transform_spfs(spfs, C)
    h = transform_onebody(h, C)
    u = transform_twobody(u, C)
    grid = copy(grid)
    
    transform = transform * C
    return SpatialSystem{typeof(basis)}(n, l, h, u, spfs, grid, basis, transform, V)
end

function System(system::PairingSystem, C)
    (; n, l, h, u) = system
    h = transform_onebody(h, C)
    u = transform_twobody(u, C)
    
    return PairingSystem(n, l, h, u)
end

function transform_spfs(spfs, C)
    spfs2 = zero.(spfs)
    l = length(spfs)
    
    for i in 1:l
        psi1 = spfs2[i]
        for j in 1:l
            psi1 .+= C[j, i] .* spfs[j]
        end
    end
    return spfs2
end

function transform_onebody(h, C)
    return C' * h * C
end

# using TensorOperations: @tensor
function transform_twobody(u, C)
    @tensor begin
        u2[a, b, c, i] := u[a, b, c, d] * C[d, i]
        u3[a, b, i, d] := u2[a, b, c, d] * C[c, i]
        u2[a, i, c, d] = u3[a, b, c, d] * C[b, i]
        u3[i, b, c, d] = u2[a, b, c, d] * C[a, i]
    end
    return u3
end