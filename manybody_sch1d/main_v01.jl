

# Kinetic operator (finite differences)
function kinetic_operator(N, dx, stencil)
    if stencil == 3
        coeffs = [1, -2, 1]
        off = (-1,0,1)
    elseif stencil == 5
        coeffs = [-1, 16, -30, 16, -1] ./ 12
        off = (-2,-1,0,1,2)
    elseif stencil == 7
        coeffs = [2, -27, 270, -490, 270, -27, 2] ./ 180
        off = (-3,-2,-1,0,1,2,3)
    elseif stencil == 9
        coeffs = [-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9] ./ 5040
        off = (-4,-3,-2,-1,0,1,2,3,4)
    elseif stencil == 11
        coeffs = [8, -125, 1000, -6000, 42000, -73766, 42000, -6000, 1000, -125, 8] ./ 25200
        off = (-5,-4,-3,-2,-1,0,1,2,3,4,5)
    elseif stencil == 13
        coeffs = [-50, 864, -7425, 44000, -222750, 1425600, -2480478, 1425600, -222750, 44000, -7425, 864, -50] ./ 831600
        off = (-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)
    else
        error("stencil must be 3,5,7,9,11,13")
    end
    s = length(coeffs) ÷ 2
    K = spdiagm(0 => fill(coeffs[s+1], N))
    for i in 1:s
        K += spdiagm(i => fill(coeffs[s+1+i], N-i))
        K += spdiagm(-i => fill(coeffs[s+1+i], N-i))
    end
    return -0.5 * K / dx^2
end

# Helper to build Kronecker product with the correct ordering:
# We want the linear index to be: idx = 1 + Σ (i_p - 1) * N^(p-1)   (particle 1 fastest)
# This corresponds to building the product as: factor N_el ⊗ ... ⊗ factor 1
function build_term(H1, I_mat, N_el, p)
    # Start with factor for particle 1 (fastest)
    term = (p == 1) ? H1 : I_mat
    for part in 2:N_el
        term = kron((part == p) ? H1 : I_mat, term)
    end
    return term
end


function solve_manybody(L, N, ω, a_soft, electrons; stencil=13, nev=10, k=0)
    # ----- 1. Grid and single‑particle operators -----
    x = range(-L, L, length=N)
    dx = step(x)
    N_el = length(electrons)

    H1 = kinetic_operator(N, dx, stencil) + spdiagm(0 => 0.5 * ω^2 * x.^2)

    # Interaction matrix (soft‑Coulomb)
    Vee = [1 / sqrt((x[i]-x[j])^2 + a_soft^2) for i in 1:N, j in 1:N]

    # ----- 2. Build many‑body Hamiltonian in the product basis -----
    dim = N^N_el
    println("Hamiltonian dimension: $dim")

    I_mat = sparse(I, N, N)
    H0 = spzeros(dim, dim)
    for p in 1:N_el
        H0 += build_term(H1, I_mat, N_el, p)
    end

    # Diagonal interaction U(x1,...,xN) = ∑_{i<j} Vee(xi, xj)
    spatial_dims = ntuple(_ -> N, N_el)
    U_diag = zeros(dim)
    for I in CartesianIndices(spatial_dims)
        idx = LinearIndices(spatial_dims)[I]   # column‑major order
        coords = Tuple(I)
        sum_val = 0.0
        for a in 1:N_el, b in a+1:N_el
            sum_val += Vee[coords[a], coords[b]]
        end
        U_diag[idx] = sum_val
    end
    H = H0 + spdiagm(0 => U_diag)

    # ----- 3. Solve for lowest eigenvalues/vectors -----
    print("Calling eigs ...")
    λ, V = eigs(H, nev=nev, which=:SR, maxiter=2000, tol=1e-10)
    println("... done")
    #
    perm = sortperm(λ)
    λ = λ[perm]
    V = V[:, perm]

    # ----- 4. Build the spin tensor (outer product of spinors) -----
    u = [1.0, 0.0]
    d = [0.0, 1.0]
    spinors = [c == 'u' ? u : d for c in electrons]
    spin = spinors[1]
    for i in 2:N_el
        spin = kron(spin, spinors[i])
    end
    spin = reshape(spin, ntuple(_ -> 2, N_el))   # shape (2,2,...,2)

    # ----- 5. Antisymmetrize each eigenstate -----
    antisym_states = []
    antisym_energies = []
    for kk in 1:nev
        #
        println("\nAntisymmetrizing state $kk")
        #
        ψ = reshape(V[:, kk], spatial_dims)      # shape (N,N,...,N)

        # Combine spatial and spin: full[x1,σ1, x2,σ2, ...]
        full_dims = ntuple(i -> isodd(i) ? N : 2, 2N_el)
        full = zeros(full_dims...)

        # Fill by iterating over all spatial and spin indices
        for Ispat in CartesianIndices(spatial_dims)
            ψ_val = ψ[Ispat]
            for Ispin in CartesianIndices(ntuple(_ -> 2, N_el))
                χ_val = spin[Ispin]
                # Build tuple of indices for full array: (x1,σ1,x2,σ2,...)
                full_idx = ntuple(2N_el) do i
                    p = (i + 1) ÷ 2
                    if isodd(i)
                        return Ispat[p]
                    else
                        return Ispin[p]
                    end
                end
                full[full_idx...] = ψ_val * χ_val
            end
        end

        # Antisymmetrize: sum over all permutations of particles with sign
        full_antisym = copy(full)
        for perm in permutations(1:N_el)
            println("perm = ", perm)
            if perm == collect(1:N_el)
                println("Skipping perm = ", perm)
                continue
            end
            # Parity: number of inversions
            parity = iseven(sum(1 for i in 1:N_el for j in i+1:N_el if perm[i] > perm[j])) ? 1 : -1
            println("parity = ", parity)
            # Build new dimension order: (x_{perm[1]}, σ_{perm[1]}, x_{perm[2]}, σ_{perm[2]}, ...)
            new_order = Int[]
            for i in 1:N_el
                push!(new_order, 2*(perm[i]-1)+1)
                push!(new_order, 2*(perm[i]-1)+2)
            end
            println("new_order = ", new_order)
            full_perm = permutedims(full, Tuple(new_order))
            full_antisym .+= parity .* full_perm
        end

        # Normalize
        norm2 = sum(abs2, full_antisym) * dx^N_el
        full_antisym ./= sqrt(norm2)

        if !isapprox(norm(full_antisym), 0, atol=1e-10)
            push!(antisym_states, full_antisym)
            push!(antisym_energies, λ[kk])
        end
    end

    println("\nAfter antisymmetrizing: antisym_energies = ", antisym_energies)

    # Remove duplicate states (same up to global phase)
    unique_states = []
    unique_energies = []
    for i in 1:length(antisym_states)
        if i == 1
            push!(unique_states, antisym_states[i])
            push!(unique_energies, antisym_energies[i])
        else
            is_new = true
            for j in 1:length(unique_states)
                if isapprox(antisym_states[i], unique_states[j], atol=1e-6) ||
                   isapprox(antisym_states[i], -unique_states[j], atol=1e-6)
                    is_new = false
                    break
                end
            end
            if is_new
                push!(unique_states, antisym_states[i])
                push!(unique_energies, antisym_energies[i])
            end
        end
    end

    if k+1 > length(unique_states)
        @warn "Requested state index $k not found, returning last available."
        idx = length(unique_states)
    else
        idx = k+1
    end

    return (ψ_antisym = unique_states[idx],
            energy = unique_energies[idx],
            all_energies = unique_energies,
            x = x,
            spin_shape = ntuple(_ -> 2, N_el),
            spatial_shape = spatial_dims)
end

# ===== Example usage =====
function main_debug()
    println("Two electrons (ud):")
    L = 5.0
    N = 60
    ω = 1.0
    a_soft = 1.0
    result = solve_manybody(L, N, ω, a_soft, "ud"; stencil=13, nev=10)
    println("Ground state energy: ", result.energy)
    println("First few energies: ", result.all_energies[1:min(5,length(result.all_energies))])

    @infiltrate

    #println("\nThree electrons (uud):")
    #N = 40
    #result3 = solve_manybody(L, N, ω, a_soft, "uud"; stencil=13, nev=10)
    #println("Ground state energy: ", result3.energy)
    #println("First few energies: ", result3.all_energies[1:min(5,length(result3.all_energies))])
end
