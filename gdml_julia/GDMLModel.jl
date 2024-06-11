mutable struct GDMLModel
    Natoms::Int64
    Ntrain::Int64
    desc_dim::Int64
    R_desc_v::Vector{Vector{Float64}}
    R_d_desc_v::Vector{Matrix{Float64}}
    σ::Float64
    λ::Float64
    R_d_desc_α::Matrix{Float64}
    indices::TrilIndices
    c::Float64
    y_std::Float64
end
