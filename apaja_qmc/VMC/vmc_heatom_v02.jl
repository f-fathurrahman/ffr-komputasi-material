using StaticArrays
import LinearAlgebra: norm
using Printf
import Random
using Statistics

include("my_vmc_stats.jl")

const D = 0.5    # hbar^2/2m in a.u. 
const Z = 2                  # charge number
const Eexact = -2.903724377  # exact ground state energy (non-relativistic, fixed nucleus)

#
# Parameters
#
const blocksize = 1_000_000   # data block size
const Ntherm = 100         # thermalization steps
const accuracy_goal = 2e-4

mutable struct VMC_Params
    Ntry::Int64
    Naccept::Int64
    step::Float64
end

mutable struct Walker
    R::Vector{MVector{3,Float64}}
    psi_ansatz::Float64
    E::Float64
end

function rand_pos(Nelectrons::Int64)
    R = Vector{MVector{3,Float64}}(undef, Nelectrons)
    for iel in 1:Nelectrons
        R[iel] = @MVector rand(3)
    end
    return R
end

# initialization 
function init()
    Nelectrons = 2   
    R = rand_pos(Nelectrons)  # coordinates
    walker = Walker(R, 0.0, 0.0) # R, psi_ansatz=0, E=0    
    vmc_params = VMC_Params(0, 0, 3.1)
    return walker, vmc_params
end

function get_wave_function_params(trial::Symbol)
    println("trial wave function params = ", trial)
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

function dist2(R::Vector{MVector{3,Float64}}, i::Int64, j::Int64)
    return norm(R[i] - R[j])
end

function dist1(R::Vector{MVector{3,Float64}}, i::Int64)
    return norm(R[i])
end

function psi_ansatz(R::Vector{MVector{3,Float64}}, wf_params)
    α, α12, β = wf_params
    r12 = dist2(R, 1, 2)
    r1 = dist1(R, 1)
    r2 = dist1(R, 2)
    b = 1 + β*r12
    return exp(-α*(r1+r2) + α12*r12/b)
end

function EL(R::Vector{MVector{3,Float64}}, wf_params)
    α, α12, β = wf_params
    r12 = dist2(R, 1, 2)
    r1  = dist1(R, 1)
    r2  = dist1(R, 2)
    b = 1 + β*r12
    ∇S = 0.5*drift(R, wf_params)
    ∇2S = -α*2/r1 -α*2/r2  + 4*α12/b^3*1/r12
    Nel = size(R, 1)
    ss = 0.0
    for i in 1:Nel
        ss += ∇S[i][1]^2 + ∇S[i][2]^2 + ∇S[i][3]^2
    end
    EL = -0.5*(ss + ∇2S) - Z/r1 - Z/r2 + 1/r12
    return EL
end

@inline function get_unitvecs( R::Vector{MVector{3,Float64}} )
    vecr1 = R[1]
    vecr2 = R[2]
    vecr12 = R[1] - R[2]
    r12 = dist2(R, 1, 2)
    r1  = dist1(R, 1)
    r2  = dist1(R, 2)
    return vecr1/r1, vecr2/r2, vecr12/r12
end

function drift(R::Vector{MVector{3,Float64}}, wf_params)
    α, α12, β = wf_params
    hatr1, hatr2, hatr12 = get_unitvecs(R)
    r12 = dist2(R, 1, 2)
    b = 1 + β*r12    
    ∇S = [-α*hatr1 + α12/b^2*hatr12, -α*hatr2 - α12/b^2*hatr12]
    return 2*∇S
end

function drift_prealloc!(R::Vector{MVector{3,Float64}}, res, wf_params)
    α, α12, β = wf_params
    hatr1, hatr2, hatr12 = get_unitvecs(R)
    r12 = dist2(R, 1, 2)
    b = 1 + β*r12
    res[1] = 2*(-α*hatr1 + α12/b^2*hatr12)
    res[2] = 2*(-α*hatr2 - α12/b^2*hatr12)
    return
end


#const buf = Vector{Float64}(undef, 1024*3) # max_N * max_dim
#const d_buf = Vector{Float64}(undef, 3)

function vmc_step!(
    R::Vector{MVector{3,Float64}}, params::VMC_Params, this_psi_ansatz::Function,
    rrs::Vector{MVector{3,Float64}}, d::MVector{3,Float64}
)
    Nelectrons = size(R, 1)
    if params.step < 1e-15
        error("VMC step is zero")
    end
    ΨR = this_psi_ansatz(R)
    if ΨR < 0.0
        error("psi_ansatz(R)<0, this shouldn't happen.")
    end
    Ψ2 = ΨR^2
    
    Random.rand!(rrs)
    #println("\nrrs in vmc_step")
    #for iel in 1:Nelectrons
    #    println(rrs[iel])
    #end
    #println("params.step = ", params.step)
    for i in 1:Nelectrons
        d[1] = params.step * ( rrs[i][1] - 0.5 )
        d[2] = params.step * ( rrs[i][2] - 0.5 )
        d[3] = params.step * ( rrs[i][3] - 0.5 )
        #println("d = ", d)
        #
        R[i][1] += d[1]
        R[i][2] += d[2]
        R[i][3] += d[3]
        Ψ_new = this_psi_ansatz(R)
        Ψ2_new = Ψ_new^2
        ratio = Ψ2_new/Ψ2
        # Metropolis:
        accept = false        
        params.Ntry += 1
        if Ψ_new < 0.0
            # fermi nodal cell change, and shoudn't happen for boson ground state either
            accept = false
        else        
            if ratio >= 1.0
                accept = true            
            else
                if ratio > rand() 
                    accept = true
                end
            end
        end
        #@info "ratio = $(ratio)"
        #@info "accept = $(accept)"
        if accept
            params.Naccept += 1
            Ψ2 = Ψ2_new     
        else
            R[i][1] -= d[1]
            R[i][2] -= d[2]
            R[i][3] -= d[3]
        end        
    end
    return sqrt(Ψ2)
end

# adjust step to keep acceptance 50-60 %
function adjust_step!(params::VMC_Params)
    minstep = 1e-5
    maxstep = 20.0
    acceptance = params.Naccept*100.0/params.Ntry
    if acceptance < 50.0
        params.step *= 0.9
    end
    if acceptance>60.0
        params.step *= 1.1
    end
    params.step = max(minstep, params.step)
    params.step = min(maxstep, params.step)
    if params.step == maxstep
        println("BAD Error: step is maxstep ",maxstep)
        @show(acceptance)
        exit()
    end
end

function eval_V(R::Vector{MVector{3,Float64}})
    r12 = dist2(R, 1, 2)
    r1  = dist1(R, 1)
    r2  = dist1(R, 2)
    return -Z/r1 - Z/r2 + 1/r12 
end


# EL check
function num_check_EL(
    R::Vector{MVector{3,Float64}}, EL_in, this_psi_ansatz::Function, this_V::Function
)
    # assumes units are D = hbar^2/(2m) = 1/2
    Ψ_in = this_psi_ansatz(R) # needed for fermion sign checks
    Rloc = copy(R)
    dim = size(R[1], 1)
    Nelectrons = size(R,1)
    h = 1e-4
    TL_num = 0.0
    psi = this_psi_ansatz(Rloc)
    for i in 1:Nelectrons
        for k in 1:dim
            Rloc[i][k] += h
            psi_p = this_psi_ansatz(Rloc)
            if psi_p*Ψ_in < 0
                Rloc = copy(R)
                println("sign change, skip")
                break
            end
            
            Rloc[i][k] -= 2h
            psi_m = this_psi_ansatz(Rloc)
            if psi_m*Ψ_in < 0
                Rloc = copy(R)
                println("sign change, skip")
                break
            end
            Rloc[i][k] += h
            TL_num += (psi_p - 2*psi + psi_m)
        end
    end
    TL_num = -1/2*TL_num/h^2 /psi
    EL_num =  TL_num + this_V(R)
    diff = (EL_in-EL_num)/EL_num
    @printf("EL check: analytical: %20.12f numerical: %20.12f  rel. difference: %20.12e\n", EL_in,EL_num,diff)
    if abs(diff)>1e-4 || isnan(diff)
        @printf("EL check: analytical: %20.12f numerical: %20.12f  rel. difference: %20.12e\n", EL_in,EL_num,diff)
        error("EL check failed")
    end
end

function num_check_∇S(
    R::Vector{MVector{3,Float64}}, this_psi_ansatz::Function, this_drift::Function
)
    S(x) = log(abs(this_psi_ansatz(x))) # psi_ansatz=exp(S)
    #
    Ψ_in = this_psi_ansatz(R) # needed for fermion sign checks
    ∇S  = 0.5*this_drift(R)  # drift = 2∇S
    dim = size(R[1], 1)
    N = size(R, 1)
    h = 1e-6
    Rloc = copy(R)
    
    for i in 1:N
        for k in 1:dim
            Rloc[i][k] += h            
            # check it's allowed, fermion psi_ansatz may change sign.
            if this_psi_ansatz(Rloc)*Ψ_in < 0
                Rloc = copy(R)
                println("sign change, skip")
                break
            end
            S_p = S(Rloc)
            Rloc[i][k] -= 2h
            if this_psi_ansatz(Rloc)*Ψ_in < 0
                Rloc = copy(R)
                println("sign change, skip")
                break 
            end
            S_m = S(Rloc)
            Rloc[i][k] += h # move back
            ∇S_num = (S_p - S_m)/(2h)
            diff = (∇S[i][k] - ∇S_num)/∇S[i][k]
            #@printf("∇S check: analytical: %20.12f numerical: %20.12f  rel difference: %20.12e\n",∇S[k,i],∇S_num,diff)
            if abs(diff) > 1e-4 || isnan(diff)
                println(k," ",i)
                @printf("∇S check: analytical: %20.12f numerical: %20.12f  rel difference: %20.12e\n",∇S[k,i],∇S_num,diff)
                error("∇S check failed")
            end
        end
    end
end

# Main program 
function main()

    Random.seed!(1234)
    Nelectrons = 2
    walker, vmc_params = init()

    # choose trial wave function parameters defined in Model_Heatom.jl
    trial = :energy_optimized_parameters
    par = get_wave_function_params(trial)
    wf_params = (par.α, par.α12, par.β) # tuple
    @show wf_params 

    # closures with fixed wf_params
    psi_ansatz_par(x) = psi_ansatz(x, wf_params)
    EL_par(x) = EL(x, wf_params)
    drift_par(x) = drift(x, wf_params)

    rrs = Vector{MVector{3,Float64}}(undef, Nelectrons)
    for i in 1:2
        rrs[i] = @MVector zeros(Float64, 3)
    end
    d_buf = @MVector zeros(Float64, 3)

    # thermalization
    println("thermalizing") 
    for i in 1:Ntherm
        vmc_step!(walker.R, vmc_params, psi_ansatz_par, rrs, d_buf)
        adjust_step!(vmc_params)
    end    
    println("thermalization done")
    
#=
    println("-"^20," extra checks ","-"^20)
    println("Checking numerically local energy EL against trial wave function psi_ansatz")
    println("to spot errors in derivatives")
    # Checks could be done more accurately using AD
    for i in 1:5
        vmc_step!(walker.R, vmc_params, psi_ansatz_par, rrs, d_buf)
        walker.E = EL_par(walker.R)
        num_check_EL(walker.R, walker.E, psi_ansatz_par, eval_V)
    end
    num_check_∇S(walker.R, psi_ansatz_par, drift_par)
    println("checks passed")
    println("-"^55)
=#

    #
    # VMC
    #
    # init E measurement (QMC_Statistics)
    Estat = init_stat(1, blocksize)
    #
    ivmc = 0
    while true # until accuracy_goal is reached
    #for iter_vmc in 1:10
        #println("Pass here iter_vmc = ", iter_vmc)
        vmc_step!(walker.R, vmc_params, psi_ansatz_par, rrs, d_buf)
        walker.E = EL_par(walker.R)

        # add new energy data
        add_sample!(Estat, walker.E)

        ivmc += 1
        if ivmc%10 == 0                      
            adjust_step!(vmc_params)
        end
        
        #
        # output when a block is full
        if Estat.finished            
            E_ave, E_std, E_inputvar2, Nb = get_stats(Estat)
            @printf("VMC %15d E = %.10f <E> = %.10f ± %.10f\n", ivmc, walker.E,  E_ave, E_std)

            if Nb > 10 && E_std < accuracy_goal
                println("reached accuracy goal")
                println("used trial wf parameters from set ",trial)
                println("Trial wf parameters:")
                @show par.α
                @show par.α12
                @show par.β
                @printf("input σ^2 = %.6f\n", E_inputvar2)
                println("result:")
                output_MCresult(E_ave, E_std)                
                break
            end                        
        end
    end

    return

end


#main()

