function integ_trapz(f_vals, grid)
    # N should be odd?
    val = sum(f_vals)
    val = val - 0.5 * (f_vals[1] + f_vals[end])
    return val * (grid[2] - grid[1])
end

function calc_onebody_integrals(basis::SpinBasis, xgrid)
    h = calc_onebody_integrals(basis.base, xgrid)
    h = kron(h, [1 0; 0 1])
    return h
end

# This is only for HarmonicOscillatorBasis
function calc_onebody_integrals(ho::HarmonicOscillatorBasis, xgrid)
    l = ho.l
    ω = ho.ω
    return Diagonal([(n + 0.5) * ω for n in 0:l-1])
end
# xgrid is not used?



# Two-body integrals for SpinBasis

function calc_twobody_integrals(basis::SpinBasis, grid, V::Interaction)
    u = calc_twobody_integrals(basis.base, grid, V)
    return _add_spin_u(u)
end


#=
When we make a spin-up and spin-down duplicate of each basis function, we get many new integrals to compute
between these new basis functions. However, most of these integrals will be zero due to opposite spins between
the basis functions, or they will be the same as the old integrals, if the spins align.
    
We loop over the latter two indices in our two-body integral, ν and λ. The matrix of elements at the indeces
u_new[:, :, ν, λ] then correspond to integrals where the first two basis functions in the integral have spins
alternating up and down.
    
'up-up'   'up-down'   up-up   ...
'down-up' 'down-down' down-up ...
 up-up     up-down    up-up   ...
 ...       ...        ...     ...
    
This 2x2 pattern of spins (that we have marked with '') repeats, and all elements in the 2x2 pattern has the same
set of spatial functions. But only one element will be non-zero, the one that has the same spins as the basis
functions numbered ν and λ.
    
The code below finds which one of these elements in the 2x2 pattern will be non-zero, and then uses this to turn
the two-body-integrals without spin into the two-body-integrals with spin.
=#
function _add_spin_u(u_old)
    l = size(u_old)[1] * 2
    u_new = zeros((l, l, l, l))
    
    for ν in 1:l
        for λ in 1:l
            if (ν % 2 == 1)     # --- UP ---
                if (λ % 2 == 1) # UP - UP
                    pattern = [1 0; 0 0]
                else            # UP - DOWN
                    pattern = [0 1; 0 0]
                end
            else                # --- DOWN ---
                if (λ % 2 == 1) # DOWN - UP
                    pattern = [0 0; 1 0]
                else            # DOWN - DOWN
                    pattern = [0 0; 0 1]
                end
            end
            
            ν_old = Int( ceil(ν / 2) )
            λ_old = Int( ceil(λ / 2) )
            @views kron!(u_new[:, :, ν, λ], u_old[:, :, ν_old, λ_old], pattern)
        end
    end
    return u_new
end



function calc_twobody_integrals(basis::SpatialBasis, xgrid, V::Interaction)
    spfs = evaluate_on_grid(basis, xgrid)
    inner = _twobody_inner_ints(spfs, xgrid, V)
    u = _twobody_outer_ints(spfs, xgrid, inner)
    return u
end

function _twobody_inner_ints(spfs, xgrid, V::Interaction)
    Nbasis = length(spfs)
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
                res = integ_trapz(fs, xgrid)
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
    #
    for κ in 1:Nbasis     
        for λ in 1:Nbasis
            @views iis[:] .= inner_ints[κ, λ, :]
            for μ in 1:Nbasis
                for ν in 1:Nbasis
                    fs .= conj.(spfs[μ]) .* iis .* spfs[ν]
                    outer_int[μ, κ, ν, λ] = integ_trapz(fs, xgrid)
                end
            end
        end
    end
    return outer_int
end

