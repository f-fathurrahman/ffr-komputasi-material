struct CCSDState{T, M}
    system::T
    mixer::M
    
    ϵ::Vector{Float64}
    f::Matrix{Float64}
    
    t1::Matrix{Float64}
    Δt1::Matrix{Float64}
    err_1::Matrix{Float64}
    
    t2::Array{Float64, 4}
    Δt2::Array{Float64, 4}
    err_2::Array{Float64, 4}
end

function CCSD(system)
    (; l, n) = system
    mixer = DIIS(l*l*n*n + l*n)
    CCSD(system, mixer)
end

function CCSD(system, mixer)
    (; l, h, u, n) = system
    
    ϵ = sp_energies(system)
    f = fock_matrix(system)
    
    Δt1 = zeros((l, n))
    t1 = zeros((l, n))
    err_1 = zeros((l, n))
    @inbounds for i in 1:n
        for a in n+1:l
            t1[a, i] = f[a, i] / (ϵ[i] - ϵ[a])
        end
    end
    
    Δt2 = zeros((l, l, n, n))
    t2 = zeros((l, l, n, n))
    err_2 = zeros((l, l, n, n))
    @inbounds for i in 1:n
        for j in 1:n 
            for a in n+1:l
                for b in n+1:l
                    t2[a, b, i, j] = u[a, b, i, j] / (ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b])
                end
            end
        end
    end
    
    return CCSDState{typeof(system), typeof(mixer)}(system, mixer, ϵ, f, t1, Δt1, err_1, t2, Δt2, err_2)
end

function energy(state::CCSDState)
    (; system) = state
    E = reference_energy(system)
    E += corr_energy(state)
    return E
end

function corr_energy(state::CCSDState)
    #=
    Returns E_CCSD - E_0
    
    where E_0 is the energy of the reference determinant
    =#
    
    (; system, f, t1, t2) = state
    (; n, l, u) = system
    
    E = 0.0
    
    @inbounds for i in 1:n
        for a in n+1:l
            E += f[i, a] * t1[a, i]
        end
    end
    
    @inbounds for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    E += 0.25 * u[i, j, a, b] * t2[a, b, i, j]
                    E += 0.5 * u[i, j, a, b] * t1[a, i] * t1[b, j]
                end
            end
        end
    end
    
    return E
end

function update!(state::CCSDState)
    (; system, mixer, ϵ, f, t1, Δt1, t2, Δt2) = state
    (; n, u) = system
    
    L = system.l # There is a collision of notation. In here, L is the number of basis functions
    
    """
    The t1 amplitude equations, from page 75 of Crawford & Schaefer
    """
    err_1 = zero(t1)
    err_2 = zero(t2)
    #@inbounds Threads.@threads 
    for a in n+1:L
        for i in 1:n
            Ω_1 = 0.0

            Ω_1 += f[a, i]
            
            for c in n+1:L
                Ω_1 += f[a, c] * t1[c, i]
            end
            
            for k in 1:n
                Ω_1 -= f[k, i] * t1[a, k]
            end

            for k in 1:n
                for c in n+1:L
                    Ω_1 += u[k, a, c, i] * t1[c, k]

                    Ω_1 += f[k, c] * t2[a, c, i, k]

                    Ω_1 -= f[k, c] * t1[c, i] * t1[a, k]
                end
            end
            øy = 0
            te = 0
            for k in 1:n
                for c in n+1:L
                    for d in n+1:L
                        Ω_1 += 0.5 * u[k, a, c, d] * t2[c, d, k, i]
                        #Ω_1 -= u[k, a, c, d] * t1[c, k] * t1[d, i] # This is from C&S, but does not work
                        Ω_1 += u[k, a, c, d] * t1[c, k] * t1[d, i] # Sign error corrected
                    end
                end
            end

            for k in 1:n
                for l in 1:n
                    for c in n+1:L
                        Ω_1 -= 0.5 * u[k, l, c, i] * t2[c, a, k, l]

                        Ω_1 -= u[k, l, c, i] * t1[c, k] * t1[a, l]
                    end
                end
            end

            for k in 1:n
                for l in 1:n
                    for c in n+1:L
                        for d in n+1:L
                            Ω_1 -= u[k, l, c, d] * t1[c, k] * t1[d, i] * t1[a, l]

                            Ω_1 += u[k, l, c, d] * t1[c, k] * t2[d, a, l, i]

                            Ω_1 -= 0.5 * u[k, l, c, d] * t2[c, d, k, i] * t1[a, l]

                            Ω_1 -= 0.5 * u[k, l, c, d] * t2[c, a, k, l] * t1[d, i]
                        end
                    end
                end
            end
            
            ϵ_ai = ϵ[i] - ϵ[a]
            Δt1[a, i] = Ω_1 / ϵ_ai
            err_1[a, i] = Ω_1
        end
    end

    """
    The t2 amplitude equations, from page 76 of Crawford & Schaefer
    """
    #@inbounds Threads.@threads 
    for a in n+1:L
        for i in 1:n
            for j in i+1:n
                for b in a+1:L
                    Ω_2 = 0.0

                    Ω_2 += u[a, b, i, j]
                    
                    for c in n+1:L
                        Ω_2 += f[b, c] * t2[a, c, i, j] # 1
                        Ω_2 -= f[a, c] * t2[b, c, i, j] # -Pab
                        
                        Ω_2 += u[a, b, c, j] * t1[c, i] # 1
                        Ω_2 -= u[a, b, c, i] * t1[c, j] # -Pij
                    end
                    
                    for k in 1:n
                        Ω_2 -= f[k, j] * t2[a, b, i, k] # -1
                        Ω_2 += f[k, i] * t2[a, b, j, k] # Pij
                        
                        Ω_2 -= u[k, b, i, j] * t1[a, k] # -1
                        Ω_2 += u[k, a, i, j] * t1[b, k] # Pab
                    end
                    
                    for k in 1:n
                        for l in 1:n
                            Ω_2 += 0.5 * u[k, l, i, j] * t2[a, b, k, l]
                            
                            Ω_2 += 0.5 * u[k, l, i, j] * t1[a, k] * t1[b, l] # 1
                            Ω_2 -= 0.5 * u[k, l, i, j] * t1[b, k] * t1[a, l] # -Pab
                        end
                    end
                    
                    for c in n+1:L
                        for d in n+1:L
                            Ω_2 += 0.5 * u[a, b, c, d] * t2[c, d, i, j]
                            
                            Ω_2 += 0.5 * u[a, b, c, d] * t1[c, i] * t1[d, j] # 1
                            Ω_2 -= 0.5 * u[a, b, c, d] * t1[c, j] * t1[d, i] # -Pij
                        end
                    end
                    
                    for k in 1:n
                        for c in n+1:L
                            Ω_2 += u[k, b, c, j] * t2[a, c, i, k] # 1
                            Ω_2 -= u[k, b, c, i] * t2[a, c, j, k] # -Pij
                            Ω_2 -= u[k, a, c, j] * t2[b, c, i, k] # -Pab
                            Ω_2 += u[k, a, c, i] * t2[b, c, j, k] # Pij Pab
                            
                            Ω_2 -= u[k, b, i, c] * t1[a, k] * t1[c, j] # -1
                            Ω_2 += u[k, b, j, c] * t1[a, k] * t1[c, i] # Pij
                            Ω_2 += u[k, a, i, c] * t1[b, k] * t1[c, j] # Pab
                            Ω_2 -= u[k, a, j, c] * t1[b, k] * t1[c, i] # -Pij Pab
                            
                            Ω_2 += f[k, c] * t1[a, k] * t2[b, c, i, j] # 1
                            Ω_2 -= f[k, c] * t1[b, k] * t2[a, c, i, j] # -Pab
                            
                            Ω_2 += f[k, c] * t1[c, i] * t2[a, b, j, k] # 1
                            Ω_2 -= f[k, c] * t1[c, j] * t2[a, b, i, k] # -Pij
                        end
                    end
                    
                    for k in 1:n
                        for l in 1:n
                            for c in n+1:L
                                Ω_2 -= u[k, l, c, i] * t1[c, k] * t2[a, b, l, j] # -1
                                Ω_2 += u[k, l, c, j] * t1[c, k] * t2[a, b, l, i] # Pij
                                
                                Ω_2 += u[k, l, i, c] * t1[a, l] * t2[b, c, j, k] # 1
                                Ω_2 -= u[k, l, j, c] * t1[a, l] * t2[b, c, i, k] # -Pij
                                Ω_2 -= u[k, l, i, c] * t1[b, l] * t2[a, c, j, k] # -Pab
                                Ω_2 += u[k, l, j, c] * t1[b, l] * t2[a, c, i, k] # Pij Pab
                                
                                Ω_2 += 0.5 * u[k, l, c, j] * t1[c, i] * t2[a, b, k, l] # 1
                                Ω_2 -= 0.5 * u[k, l, c, i] * t1[c, j] * t2[a, b, k, l] # -Pij
                                
                                Ω_2 += 0.5 * u[k, l, c, j] * t1[c, i] * t1[a, k] * t1[b, l] # 1
                                Ω_2 -= 0.5 * u[k, l, c, i] * t1[c, j] * t1[a, k] * t1[b, l] # -Pij
                                Ω_2 -= 0.5 * u[k, l, c, j] * t1[c, i] * t1[b, k] * t1[a, l] # -Pab
                                Ω_2 += 0.5 * u[k, l, c, i] * t1[c, j] * t1[b, k] * t1[a, l] # Pij Pab
                            end
                        end
                    end
                    
                    for k in 1:n
                        for c in n+1:L
                            for d in n+1:L
                                Ω_2 += u[k, a, c, d] * t1[c, k] * t2[d, b, i, j] # 1
                                Ω_2 -= u[k, b, c, d] * t1[c, k] * t2[d, a, i, j] # -Pab
                                
                                Ω_2 += u[a, k, d, c] * t1[d, i] * t2[b, c, j, k] # 1
                                Ω_2 -= u[a, k, d, c] * t1[d, j] * t2[b, c, i, k] # -Pij
                                Ω_2 -= u[b, k, d, c] * t1[d, i] * t2[a, c, j, k] # -Pab
                                Ω_2 += u[b, k, d, c] * t1[d, j] * t2[a, c, i, k] # Pij Pab
                                
                                Ω_2 -= 0.5 * u[k, b, c, d] * t1[a, k] * t2[c, d, i, j] # -1
                                Ω_2 += 0.5 * u[k, a, c, d] * t1[b, k] * t2[c, d, i, j] # Pab
                                
                                Ω_2 -= 0.5 * u[k, b, c, d] * t1[c, i] * t1[a, k] * t1[d, j] # -1
                                Ω_2 += 0.5 * u[k, b, c, d] * t1[c, j] * t1[a, k] * t1[d, i] # Pij
                                Ω_2 += 0.5 * u[k, a, c, d] * t1[c, i] * t1[b, k] * t1[d, j] # Pab
                                Ω_2 -= 0.5 * u[k, a, c, d] * t1[c, j] * t1[b, k] * t1[d, i] # -Pij Pab
                            end
                        end
                    end
                    
                    for k in 1:n
                        for l in 1:n
                            for c in n+1:L
                                for d in n+1:L
                                    Ω_2 += 0.5 * u[k, l, c, d] * t2[a, c, i, k] * t2[d, b, l, j] # 1
                                    Ω_2 -= 0.5 * u[k, l, c, d] * t2[a, c, j, k] * t2[d, b, l, i] # -Pij
                                    Ω_2 -= 0.5 * u[k, l, c, d] * t2[b, c, i, k] * t2[d, a, l, j] # -Pab
                                    Ω_2 += 0.5 * u[k, l, c, d] * t2[b, c, j, k] * t2[d, a, l, i] # Pij Pab
                                    
                                    Ω_2 += 0.25 * u[k, l, c, d] * t2[c, d, i, j] * t2[a, b, k, l]

                                    Ω_2 -= 0.5 * u[k, l, c, d] * t2[a, c, i, j] * t2[b, d, k, l] # -1
                                    Ω_2 += 0.5 * u[k, l, c, d] * t2[b, c, i, j] * t2[a, d, k, l] # Pab
                                    
                                    Ω_2 -= 0.5 * u[k, l, c, d] * t2[a, b, i, k] * t2[c, d, j, l] # -1
                                    Ω_2 += 0.5 * u[k, l, c, d] * t2[a, b, j, k] * t2[c, d, i, l] # Pij
                                    
                                    Ω_2 -= u[k, l, c, d] * t1[c, k] * t1[d, i] * t2[a, b, l, j] # -1
                                    Ω_2 += u[k, l, c, d] * t1[c, k] * t1[d, j] * t2[a, b, l, i] # Pij
                                    
                                    Ω_2 -= u[k, l, c, d] * t1[c, k] * t1[a, l] * t2[d, b, i, j] # -1
                                    Ω_2 += u[k, l, c, d] * t1[c, k] * t1[b, l] * t2[d, a, i, j] # Pab
                                    
                                    Ω_2 += 0.25 * u[k, l, c, d] * t1[c, i] * t1[d, j] * t2[a, b, k, l] # 1
                                    Ω_2 -= 0.25 * u[k, l, c, d] * t1[c, j] * t1[d, i] * t2[a, b, k, l] # -Pij
                                    
                                    Ω_2 += 0.25 * u[k, l, c, d] * t1[a, k] * t1[b, l] * t2[c, d, i, j] # 1
                                    Ω_2 -= 0.25 * u[k, l, c, d] * t1[b, k] * t1[a, l] * t2[c, d, i, j] # -Pab
                                    
                                    Ω_2 += u[k, l, c, d] * t1[c, i] * t1[b, l] * t2[a, d, k, j] # 1
                                    Ω_2 -= u[k, l, c, d] * t1[c, j] * t1[b, l] * t2[a, d, k, i] # -Pij
                                    Ω_2 -= u[k, l, c, d] * t1[c, i] * t1[a, l] * t2[b, d, k, j] # -Pab
                                    Ω_2 += u[k, l, c, d] * t1[c, j] * t1[a, l] * t2[b, d, k, i] # Pij Pab
                                    
                                    Ω_2 += 0.25 * u[k, l, c, d] * t1[c, i] * t1[a, k] * t1[d, j] * t1[b, l] # 1
                                    Ω_2 -= 0.25 * u[k, l, c, d] * t1[c, j] * t1[a, k] * t1[d, i] * t1[b, l] # -Pij
                                    Ω_2 -= 0.25 * u[k, l, c, d] * t1[c, i] * t1[b, k] * t1[d, j] * t1[a, l] # -Pab
                                    Ω_2 += 0.25 * u[k, l, c, d] * t1[c, j] * t1[b, k] * t1[d, i] * t1[a, l] # Pij Pab
                                end
                            end
                        end
                    end
                    
                    ϵ_abij = ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b]
                    Δt2_abij = Ω_2 / ϵ_abij
                    
                    Δt2[a, b, i, j] = Δt2_abij
                    Δt2[a, b, j, i] = -Δt2_abij
                    Δt2[b, a, i, j] = -Δt2_abij
                    Δt2[b, a, j, i] = Δt2_abij
                    
                    err_2[a, b, i, j] = Ω_2
                    err_2[a, b, j, i] = -Ω_2
                    err_2[b, a, i, j] = -Ω_2
                    err_2[b, a, j, i] = Ω_2 
                end
            end
        end
    end
    trial = vcat(vec(t1), vec(t2))
    direction = vcat(vec(Δt1), vec(Δt2))
    err = vcat(vec(err_1), vec(err_2))
    
    new_ts = compute_new_vector(mixer, trial, direction, err)
    t1 .= reshape(new_ts[1:L*n], size(t1))
    t2 .= reshape(new_ts[L*n+1:end], size(t2))
    
    return state
end
;