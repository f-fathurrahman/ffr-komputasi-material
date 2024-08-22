function trapz(f_vals, grid)
    # N should be odd?
    val = sum(f_vals)
    val = val - 0.5 * (f_vals[1] + f_vals[end])
    return val * (grid[2] - grid[1])
end

function calc_onebody_integrals(ho::HarmonicOscillatorBasis, xgrid)
    l = ho.l
    ω = ho.ω
    return Diagonal([(n + 0.5) * ω for n in 0:l-1])
end
# xgrid is not used?


function calc_onebody_integrals(basis::SpinBasis, xgrid)
    h = calc_onebody_integrals(basis.base, xgrid)
    h = kron(h, [1 0; 0 1])
    return h
end

function calc_twobody_integrals(spfs, xgrid, V::Interaction)
    inner = _twobody_inner_ints(spfs, xgrid, V)
    u = _twobody_outer_ints(spfs, xgrid, inner)
    return u
end

function _twobody_inner_ints(spfs, xgrid, V::Interaction)
    Nbasis = length(spfs) # including spin
    Npoints = size(spfs[1], 1)
    @assert size(xgrid,1) == Npoints
    inner_int = zeros(Float64, Nbasis, Nbasis, length(spfs[1]))
    V_tmp = zeros(Float64, Npoints)
    fs = zeros(Float64, Npoints)
    for i in 1:Npoints
        xi = xgrid[i]
        evaluate_on_grid!(V_tmp, xi, xgrid, V)
        for κ in 1:Nbasis
            for λ in κ:Nbasis
                fs .= conj.(spfs[κ]) .* V_tmp .* spfs[λ]
                res = trapz(fs, xgrid)
                inner_int[κ, λ, i] = res
                inner_int[λ, κ, i] = res'
            end
        end
    end
    return inner_int
end


function _twobody_outer_ints(spfs, xgrid, inner_ints)
    Nbasis = size(spfs, 1)
    outer_int = zeros(Float64, Nbasis, Nbasis, Nbasis, Nbasis)
    Npoints = size(xgrid, 1)
    fs = zeros(Float64, Npoints)
    iis = zeros(Float64, Npoints)

    @inbounds for κ in 1:Nbasis     
        for λ in 1:Nbasis
            @views iis[:] .= inner_ints[κ, λ, :]
            for μ in 1:Nbasis
                for ν in 1:Nbasis
                    fs .= conj.(spfs[μ]) .* iis .* spfs[ν]
                    outer_int[μ, κ, ν, λ] = trapz(fs, xgrid)
                end
            end
        end
    end
    return outer_int
end

