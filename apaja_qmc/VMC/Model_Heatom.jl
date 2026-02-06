module Model_Heatom

using StaticArrays
using LinearAlgebra: norm


push!(LOAD_PATH,".")
using Utilities: dist1, dist2



export EL, drift, ln_psi2, Ψ, N, dim, Eexact,  λ
export V
export trial, α, α12, β
export Elocal_multi, ln_psi2_multi, drift_multi, psi_multi
export Ψ_i_per_Ψ_multi
export get_wave_function_params

const λ = 0.5   # hbar^2/(2m) in a.u.

const N = 2                  # number of electrons
const Z = 2                  # charge number
const dim = 3                # dimension  
const Eexact = -2.903724377  # exact ground state energy (non-relativistic, fixed nucleus)

# =================
# choose trial type
# =================
   

function get_wave_function_params(trial::Symbol)
    @show trial
    # returns named tuples
    if trial === :energy_optimized_1S
        # energy optimized, analytical:
        # 27/16 = 1.6875
        # exact E = α^2 - 27/8*α 
        return (α = 27/16, α12 = 0.0, β = 0.0)
    elseif trial === :cusp_condition_parameters 
        # cusp condition values:
        return (α = 2.0, α12 = 0.5, β = 0.0)
    elseif trial === :three_parameters 
        return (α = 1.0, α12 = 0.5, β = 0.1)
    elseif trial === :initial_parameters_for_optimization
        return (α = 2.0, α12 = 0.5, β = 0.2)
    elseif trial === :energy_optimized_parameters
        # <E> = -2.891115 +/- 0.000020 
        return (α = 1.847529, α12 = 0.359070, β = 0.159321)
    elseif trial === :variance_optimized_parameters
        return (α = 1.949544, α12 = 0.526500, β = 0.430728)
    else
        error("Unknown trial: $trial")
    end
end


# utility
@inline function get_unitvecs(R::MMatrix{dim,N,Float64})
    vecr1 = R[:,1]
    r1 = norm(vecr1)  
    vecr2 = R[:,2]
    r2 = norm(vecr2)
    hatr2  = vecr2/r2
    vecr12 =  R[:,1]-R[:,2]
    r12 = norm(vecr12)
    vecr1/r1,vecr2/r2,vecr12/r12
end

# potential energy
@inline function V(R::MMatrix{dim,N,Float64})
    r12 = norm(R[:,1]-R[:,2])
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    - Z/r1 - Z/r2 + 1/r12 
end


#=
# Fixed wf parameter codes
# ========================
function EL(R ::MMatrix)
    r12 = norm(R[:,1]-R[:,2])
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    hatr1, hatr2, hatr12 = get_unitvecs(R)
    b = 1+β*r12
    ∇S = 0.5*drift(R)
    ∇2S = -α*2/r1 -α*2/r2  + 4*α12/b^3*1/r12 
    EL = -0.5*(sum(∇S.^2) + ∇2S) - Z/r1 - Z/r2 + 1/r12
end


# 2∇S
@inline function drift(R ::MMatrix)
    hatr1, hatr2, hatr12 =  get_unitvecs(R)
    r12 = norm(R[:,1]-R[:,2])
    b = 1+β*r12    
    ∇S = hcat(-α*hatr1 + α12/b^2*hatr12, -α*hatr2 - α12/b^2*hatr12)
    2∇S
end

# ln(φ_T^2) = 2*ln(φ_T)
function ln_psi2(R ::MMatrix)
    r12 = norm(R[:,1]-R[:,2])
    r1 = norm(R[:,1])
    r2 = norm(R[:,2])
    b = 1+β*r12
    2*(-α*(r1+r2) + α12*r12/b)
end

function Ψ(R ::MMatrix)
    r12 = norm(R[:,1]-R[:,2])
    r1 = norm(R[:,1])
    r2 = norm(R[:,2])
    b = 1+β*r12
    exp(-α*(r1+r2) + α12*r12/b)
end
=#

# three-parameter functions

function Ψ(R::MMatrix{dim, N, Float64}, wf_params) where {dim, N}
    α, α12, β = wf_params
    r12 = dist2(R, 1, 2)
    r1 = dist1(R, 1)
    r2 = dist1(R, 2)
    b = 1+β*r12
    return exp(-α*(r1+r2) + α12*r12/b)
end


function EL(R::MMatrix, wf_params)
    α, α12, β = wf_params
    r12 = norm(R[:,1]-R[:,2])
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    hatr1, hatr2, hatr12 = get_unitvecs(R)
    b = 1+β*r12
    ∇S = 0.5*drift(R, wf_params)
    ∇2S = -α*2/r1 -α*2/r2  + 4*α12/b^3*1/r12
    EL = -0.5*(sum(∇S.^2) + ∇2S) - Z/r1 - Z/r2 + 1/r12
end

# 2∇S
@inline function drift(R::MMatrix, wf_params)
    α, α12, β = wf_params
    hatr1, hatr2, hatr12 =  get_unitvecs(R)
    r12 = norm(R[:,1]-R[:,2])
    b = 1+β*r12    
    ∇S = hcat(-α*hatr1 + α12/b^2*hatr12, -α*hatr2 - α12/b^2*hatr12)
    2∇S
end


# Multiparameter functions
# ------------------------
# parameters c, α1, α2, α12, β
#
# ψ :=  sum_k ck^2*exp(-α1*r1 - α2*r2 - α12*r12)  + (1<->2)
#   :=  sum_k ck^2*exp(Sk) + (1<->2) 
#   Sk :=  -α1*r1 - α2*r2 - α12*r12 
#    
# ln(|Ψ|^2) = 2ln|Ψ| = 2*ln( sum_k ck*exp(Sk) + (1<->2) )
# ∇Ψ_i =  sum_k ck*exp(Sk) ∇_iSk +  (1<->2)  
# drift_i := 2(∇_iΨ)/Ψ = 2/Ψ* [ sum_k ck*exp(Sk) ∇Sk  + (1<->2) ]
#       
# ∇_i^2Ψ =  sum_k ck^2*exp(Sk) (∇_i^2Sk + (∇_iSk)^2)  + (1<->2) 
# TL = -1/2*Σ_i (∇_i^2Ψ)/Ψ = -1/(2Ψ)* [ sum_k ck^2*exp(Sk) Σ_i(∇_i^2Sk + (∇_iSk)^2)  + (1<->2) ] 
# EL = TL + V(R)
#
@inline function Elocal_multi(R ::MMatrix, wf_params ::Vector{Float64})    
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    r12 = norm(R[:,1]-R[:,2])
    hatr1, hatr2, hatr12 =  get_unitvecs(R)
    
    TL = 0.0
    Ψ = 0.0    
    for k in 1:4:length(wf_params)
        ck, α1, α2, α12 = wf_params[k:k+3]
        if abs(ck)<1e-15 continue end
        for sym = 1:2
            α1, α2 =  α2, α1                   
            Sk = -α1*r1 - α2*r2 - α12*r12
            ∇Sk = hcat(-α1*hatr1 - α12*hatr12, -α2*hatr2 + α12*hatr12)
            ∇2Sk = -α1*2/r1 -α2*2/r2  - 4*α12*1/r12
            Ψ  += ck^2*exp(Sk)
            TL += ck^2*exp(Sk)* (sum(∇Sk.^2) + ∇2Sk)
        end
    end
    TL *= -1/(2Ψ)
    EL = TL + V(R)
end



# drift_i := 2∇_iΨ/Ψ = 2/Ψ* [ sum_k ck*exp(Sk) ∇_iSk  + (1<->2) ]
# ∇_iΨ =  sum_k ck*exp(Sk) ∇_iSk +  (1<->2)
# Careful not to mix ∇_1 Ψ and ∇_2 Ψ, the gradients are *not* supposed to be symmetrized
@inline function drift_multi(R ::MMatrix, wf_params ::Vector{Float64})
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    r12 = norm(R[:,1]-R[:,2])
    hatr1, hatr2, hatr12 =  get_unitvecs(R)
    
    ∇Ψ = zeros(MMatrix{3,2})
    Ψ = 0.0
    for sym = 1:2
        r1, r2 = r2, r1
        for k in 1:4:length(wf_params)
            ck, α1, α2, α12 = wf_params[k:k+3]
            Sk = -α1*r1 - α2*r2 - α12*r12
            Ψ += ck^2*exp(Sk)
            if sym==2
                ∇Ψ += ck*exp(Sk)* hcat(-α1*hatr1 - α12*hatr12, -α2*hatr2 + α12*hatr12)
            else
                ∇Ψ += ck*exp(Sk)* hcat(-α2*hatr1 - α12*hatr12, -α1*hatr2 + α12*hatr12)
            end
        end
    end

    
    drift = 2∇Ψ/Ψ 
    drift
end



# ψ = [  sum_k ck*exp(-α1*r1 - α2*r2 - α12*r12)  + (1<->2) ]
#   := [ sum_k ck*exp(Sk) + (1<->2) ]
#   Sk :=  -α1*r1 - α2*r2 - α12*r12 
@inline function psi_multi(R ::MMatrix, wf_params ::Vector{Float64})
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    r12 = norm(R[:,1]-R[:,2])
    Ψ = 0.0
    for sym = 1:2
        r1, r2 = r2, r1
        for k in 1:4:length(wf_params)
            ck, α1, α2, α12 = wf_params[k:k+3]
            Sk = -α1*r1 - α2*r2 - α12*r12 
            Ψ += ck^2*exp(Sk)
        end
    end
    Ψ
end

# ψ = [  sum_k ck*exp(-α1*r1 - α2*r2 -α12*r12)  + (1<->2) ]
#   := [ sum_k ck*exp(Sk) + (1<->2) ]
#   Sk :=  -α1*r1 - α2*r2 - α12*r12 
@inline function ln_psi2_multi(R ::MMatrix, wf_params ::Vector{Float64})
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    r12 = norm(R[:,1]-R[:,2])
    npara = length(wf_params)
    Ψ = 0.0
    for sym = 1:2
        r1, r2 = r2, r1
        for k in 1:4:npara
            ck, α1, α2, α12 = wf_params[k:k+3]
            if abs(ck)<1e-15 continue end
            Sk = -α1*r1 - α2*r2 - α12*r12
            Ψ += ck^2*exp(Sk)
        end
    end   
    2*log(Ψ) # will give error if Ψ<0, but boson ground state wf is non-negative 
end

# Ψ_i/Ψ 
@inline function Ψ_i_per_Ψ_multi(R ::MMatrix, wf_params ::Vector{Float64})
    r1  = norm(R[:,1])
    r2  = norm(R[:,2])
    r12 = norm(R[:,1]-R[:,2])
    npara = length(wf_params)
    Ψ = 0.0    
    Ψi = zeros(npara)
    for sym = 1:2
        r1, r2 = r2, r1
        for k in 1:4:npara
            ck, α1, α2, α12 = wf_params[k:k+3]
            Sk = -α1*r1 - α2*r2 - α12*r12
            ee = exp(Sk)
            Ψ += ck^2*ee
            Ψi[k] += 2ck*ee
            Ψi[k+1] += -r1 * ck^2*ee
            Ψi[k+2] += -r2 * ck^2*ee
            Ψi[k+3] += -r12 * ck^2*ee
        end
    end
    Ψi/Ψ

end

end




