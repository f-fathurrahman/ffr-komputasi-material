struct CCDState{T, M}
    system::T
    mixer::M
    
    ϵ::Vector{Float64}
    f::Matrix{Float64}
    
    t::Array{Float64, 4}
    Δt::Array{Float64, 4}
    err::Array{Float64, 4}
end

function CCD(system)
    (; l, n) = system
    mixer = DIIS((l, l, n, n))
    CCD(system, mixer)
end

function CCD(system, mixer)
    (; l, u, n) = system
    
    ϵ = sp_energies(system)
    f = fock_matrix(system)
    
    t = zeros((l, l, n, n))
    Δt = zeros((l, l, n, n))
    err = zeros((l, l, n, n))
    
    @inbounds for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    t[a, b, i, j] = u[a, b, i, j] / (ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b])
                end
            end
        end
    end
    
    return CCDState{typeof(system), typeof(mixer)}(system, mixer, ϵ, f, t, Δt, err)
end

function energy(state::CCDState)
    (; system) = state
    
    E = reference_energy(system)
    E += corr_energy(state)
    return E
end

function corr_energy(state::CCDState)
    (; system, t) = state
    (; n, l, u) = system
    
    E = 0.0
    @inbounds for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    E += 0.25 * u[i, j, a, b] * t[a, b, i, j]
                end
            end
        end
    end
    return E
end

function update!(state::CCDState)
    (; system, mixer, ϵ, f, t, Δt, err) = state
    (; n, l, u) = system
    
    err .= zero(Float64)
    Threads.@threads for a in n+1:l
        for i in 1:n
            for j in i+1:n
                for b in a+1:l
                    # Page 361 of An Advanced Course in Computational Physics
                    Ω = 0
                    
                    for c in n+1:l
                        Ω += f[b, c] * t[a, c, i, j] # 1
                        Ω -= f[a, c] * t[b, c, i, j] # -Pab
                    end
                    
                    for k in 1:n
                        Ω -= f[k, j] * t[a, b, i, k] # -1
                        Ω += f[k, i] * t[a, b, j, k] # Pij
                    end
                    
                    Ω += u[a, b, i, j]
                    
                    for c in n+1:l
                        for d in n+1:l
                            Ω += 0.5 * u[a, b, c, d] * t[c, d, i, j]
                        end
                    end

                    for k in 1:n
                        for li in 1:n
                            Ω += 0.5 * u[k, li, i, j] * t[a, b, k, li]
                        end
                    end

                    for k in 1:n
                        for c in n+1:l
                            Ω += u[k, b, c, j] * t[a, c, i, k] # 1
                            Ω -= u[k, b, c, i] * t[a, c, j, k] # -Pij
                            Ω -= u[k, a, c, j] * t[b, c, i, k] # -Pab
                            Ω += u[k, a, c, i] * t[b, c, j, k] # Pij Pab
                        end
                    end

                    for k in 1:n
                        for li in 1:n
                            for c in n+1:l
                                for d in n+1:l
                                    Ω += 0.25 * u[k, li, c, d] * t[c, d, i, j] * t[a, b, k, li]

                                    Ω += u[k, li, c, d] * t[a, c, i, k] * t[b, d, j, li] # 1
                                    Ω -= u[k, li, c, d] * t[a, c, j, k] * t[b, d, i, li] # -Pij

                                    Ω -= 0.5 * u[k, li, c, d] * t[d, c, i, k] * t[a, b, li, j] # -1
                                    Ω += 0.5 * u[k, li, c, d] * t[d, c, j, k] * t[a, b, li, i] # Pij

                                    Ω -= 0.5 * u[k, li, c, d] * t[a, c, li, k] * t[d, b, i, j] # -1
                                    Ω += 0.5 * u[k, li, c, d] * t[b, c, li, k] * t[d, a, i, j] # Pab
                                end
                            end
                        end
                    end
                    
                    ϵ_abij = ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b]
                    Δt_abij = Ω / ϵ_abij

                    Δt[a, b, i, j] = Δt_abij
                    Δt[a, b, j, i] = -Δt_abij
                    Δt[b, a, i, j] = -Δt_abij
                    Δt[b, a, j, i] = Δt_abij
                    
                    err[a, b, i, j] = Ω
                    err[a, b, j, i] = -Ω
                    err[b, a, i, j] = -Ω
                    err[b, a, j, i] = Ω    
                end
            end
        end
    end
    
    t .= compute_new_vector(mixer, t, Δt, err)
    
    return state
end
;