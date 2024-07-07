abstract type Mixer end

mutable struct Alpha <: Mixer
    α::Float64
    function Alpha(α)
        @assert 0 <= α <= 1
        return new(α)
    end
end

function compute_new_vector(mixer::Alpha, t, Δt, error)
    (; α) = mixer
    return t .+ Δt .* α
end

mutable struct DIIS <: Mixer
    iter::Int
    
    max_vecs::Int
    
    trial_vecs::Vector{Vector{Float64}}
    direction_vecs::Vector{Vector{Float64}}
    error_vecs::Vector{Vector{Float64}}
end

function DIIS(shape, max_vecs = 10)
    vec_length = prod(shape)
    trial_vecs = [zeros(vec_length) for i in 1:max_vecs]
    direction_vecs = [zeros(vec_length) for i in 1:max_vecs]
    error_vecs = [zeros(vec_length) for i in 1:max_vecs]
    
    return DIIS(0, max_vecs, trial_vecs, direction_vecs, error_vecs)
end

function compute_new_vector(mixer::DIIS, trial, direction, error)
    (; iter, max_vecs, trial_vecs, direction_vecs, error_vecs) = mixer
    
    index = iter % max_vecs + 1
    iter += 1
    mixer.iter = iter
    
    error_vecs[index] .= vec(error)
    trial_vecs[index] .= vec(trial)
    direction_vecs[index] .= vec(direction)
    
    # Setting up the equations
    B_dim = min(iter, max_vecs) + 1
    B = zeros(B_dim, B_dim)
    for i in 1:B_dim-1
        for j in 1:i
            B[i, j] = LinearAlgebra.dot(error_vecs[i], error_vecs[j])
            if i != j
                B[j, i] = B[i, j]
            end
        end
        B[i, B_dim] = -1
        B[B_dim, i] = -1
    end
    
    pre_condition = zeros(B_dim)
    if any(i <= 0 for i in LinearAlgebra.diag(B)[1:end-1])
        pre_condition[1:end-1] .= 1
    else
        pre_condition[1:end-1] .= 1 ./ sqrt.(LinearAlgebra.diag(B)[1:end-1])
    end
    pre_condition[end] = 1
    
    for i in 1:B_dim
        for j in 1:B_dim
            B[i, j] *= pre_condition[i] * pre_condition[j]
        end
    end
    # Solving the equations
    weights = -LinearAlgebra.pinv(B)[end, :] # final row
    weights .*= pre_condition
    
    # Using solution to get new vector
    new_trial = zero(trial_vecs[1])
    for i in 1:B_dim-1
        new_trial .+= weights[i] .* (trial_vecs[i] .+ direction_vecs[i])
    end
    
    return reshape(new_trial, size(trial))
end