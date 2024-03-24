# Let's put all variables here for the moment
mutable struct GlobalVariables
    Natoms::Int64
    Nelectrons::Int64
    RK::Matrix{Float64}
    RE::Matrix{Float64}
    RNEU::Matrix{Float64}
end

function setup_global_vars!(global_vars)

end