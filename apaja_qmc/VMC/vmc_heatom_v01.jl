#=
JULIA

VMC on a helium atom with trial wave function φ_T
defined in module Model_Heatom

The Hamiltonian is
H = -1/2* (nabla_1^2 + nabla_2^2) - 2/r1 - 2/r2 + 1/r12  
=#

import LinearAlgebra: norm
using Printf
using Random: rand!
import Random
using Statistics

# local modules:
push!(LOAD_PATH, ".")

const dim = 3    # space dimension
const D = 0.5    # hbar^2/2m in a.u. 
const N = 2                  # number of electrons
const Z = 2                  # charge number
const Eexact = -2.903724377  # exact ground state energy (non-relativistic, fixed nucleus)

#
# Parameters
#
const blocksize = 1_000_000   # data block size
const Ntherm = 100         # thermalization steps
const accuracy_goal = 2e-4

mutable struct VMC_Params
    Ntry :: Int64
    Naccept :: Int64
    step :: Float64
end

mutable struct Walker
    R::Matrix{Float64}
    ψ::Float64
    E::Float64
end

mutable struct t_StatData
    n    ::Int64  
    data ::Vector{Float64}
    data2 ::Vector{Float64}
    input_σ2 :: Float64
end

mutable struct t_Stat    
    nblocks   ::Int64
    blocksize ::Int64
    finished  ::Bool
    sample    ::t_StatData
    datablock :: Vector{t_StatData}
    t_Stat() = new()
end

# initialization 
function init()     
    # 
    # Initialize walker
    #
    R = rand(Float64, dim, N)    # coordinates
    walker = Walker(R, 0.0, 0.0) # R, ψ=0, E=0    
    vmc_params = VMC_Params(0, 0, 3.1)
    return walker, vmc_params
end

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

function dist2(R::Matrix{Float64}, i::Int64, j::Int64)
    s = 0.0
    for k in 1:dim
        dx = R[k,i] - R[k,j]
        s += dx*dx
    end
    return sqrt(s)
end

function dist1(R::Matrix{Float64}, i::Int64)
    s = 0.0
    for k in 1:dim
        s += R[k,i]^2
    end
    return sqrt(s)
end

function Ψ(R::Matrix{Float64}, wf_params)
    α, α12, β = wf_params
    r12 = dist2(R, 1, 2)
    r1 = dist1(R, 1)
    r2 = dist1(R, 2)
    b = 1+β*r12
    return exp(-α*(r1+r2) + α12*r12/b)
end

function EL(R::Matrix{Float64}, wf_params)
    α, α12, β = wf_params
    @views r12 = norm(R[:,1]-R[:,2])
    @views r1  = norm(R[:,1])
    @views r2  = norm(R[:,2])
    b = 1 + β*r12
    ∇S = 0.5*drift(R, wf_params)
    ∇2S = -α*2/r1 -α*2/r2  + 4*α12/b^3*1/r12
    EL = -0.5*(sum(∇S.^2) + ∇2S) - Z/r1 - Z/r2 + 1/r12
    return EL
end

function get_unitvecs( R::Matrix{Float64} )
    @views vecr1 = R[:,1]
    @views vecr2 = R[:,2]
    @views vecr12 =  R[:,1] - R[:,2]
    r12 = sqrt( (R[1,1]-R[1,2])^2 + (R[2,1]-R[2,2])^2 + (R[3,1]-R[3,2])^2 )
    r1  = sqrt( R[1,1]^2 + R[2,1]^2 + R[3,1]^2 )
    r2  = sqrt( R[1,2]^2 + R[2,2]^2 + R[3,2]^2 )
    return vecr1/r1, vecr2/r2, vecr12/r12
end

function drift(R::Matrix{Float64}, wf_params)
    α, α12, β = wf_params
    hatr1, hatr2, hatr12 = get_unitvecs(R)
    @views r12 = norm(R[:,1]-R[:,2])
    b = 1+β*r12    
    ∇S = hcat(-α*hatr1 + α12/b^2*hatr12, -α*hatr2 - α12/b^2*hatr12)
    return 2*∇S
end

#const buf = Vector{Float64}(undef, 1024*3) # max_N * max_dim
#const d_buf = Vector{Float64}(undef, 3)

function vmc_step!(
    R::Matrix{Float64}, params::VMC_Params, Ψ::Function,
    rrs::Matrix{Float64}, d::Vector{Float64}
)

    #  Ψ must be callable (function or functor)
    # sanity:
    @assert N <= 1024
    if params.step < 1e-15
        error("VMC step is zero")
    end
    ΨR = Ψ(R)
    if ΨR < 0.0
        error("Ψ(R)<0, this shouldn't happen.")
    end
    Ψ2 = ΨR^2
    
    rand!(rrs)
    for i in 1:N
        #@views rr = rrs[:, i]
        #@. d  = params.step*( rr - 0.5 )
        d[1] = params.step * (rrs[1,i] - 0.5)
        d[2] = params.step * (rrs[2,i] - 0.5)
        d[3] = params.step * (rrs[3,i] - 0.5)
        #
        R[1,i] += d[1]
        R[2,i] += d[2]
        R[3,i] += d[3]
        Ψ_new = Ψ(R)
        Ψ2_new = Ψ_new^2
        ratio = Ψ2_new/Ψ2
        # Metropolis:
        accept = false        
        params.Ntry +=1
        if Ψ_new<0.0
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
        if accept
            params.Naccept += 1
            Ψ2 = Ψ2_new     
        else
            # revert to old value
            R[1,i] -= d[1]
            R[2,i] -= d[2]
            R[3,i] -= d[3]
        end        
    end
    return sqrt(Ψ2)
end

# adjust step to keep acceptance 50-60 %
function adjust_step!(params ::VMC_Params)
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

function eval_V(R::Matrix{Float64})
    #@views r12 = norm(R[:,1]-R[:,2])
    r12 = sqrt( (R[1,1]-R[1,2])^2 + (R[2,1]-R[2,2])^2 + (R[3,1]-R[3,2])^2 )
    r1  = sqrt( R[1,1]^2 + R[2,1]^2 + R[3,1]^2 )
    r2  = sqrt( R[1,2]^2 + R[2,2]^2 + R[3,2]^2 )
    return -Z/r1 - Z/r2 + 1/r12 
end


# EL check
function num_check_EL(R::Matrix{Float64}, EL_in, Ψ ::Function, V::Function)
    # assumes units are D = hbar^2/(2m) = 1/2
    Ψ_in = Ψ(R) # needed for fermion sign checks
    Rloc = copy(R)
    dim = size(R,1)
    N = size(R,2)
    h = 1e-4
    TL_num = 0.0
    psi = Ψ(Rloc)
    for i in 1:N
        for k in 1:dim
            Rloc[k,i] += h
            psi_p = Ψ(Rloc)
            if psi_p*Ψ_in<0
                Rloc = copy(R)
                println("sign change, skip")
                break
            end
            
            Rloc[k,i] -= 2h
            psi_m = Ψ(Rloc)
            if psi_m*Ψ_in<0
                Rloc = copy(R)
                println("sign change, skip")
                break
            end
            Rloc[k,i] += h
            TL_num += (psi_p - 2*psi + psi_m)
        end
    end
    TL_num = -1/2*TL_num/h^2 /psi
    EL_num =  TL_num + V(R)
    diff = (EL_in-EL_num)/EL_num
    @printf("EL check: analytical: %20.12f numerical: %20.12f  rel. difference: %20.12e\n", EL_in,EL_num,diff)
    if abs(diff)>1e-4 || isnan(diff)
        @printf("EL check: analytical: %20.12f numerical: %20.12f  rel. difference: %20.12e\n", EL_in,EL_num,diff)
        error("EL check failed")
    end
end

function num_check_∇S(R::Matrix, Ψ::Function, drift::Function)
    S(x) = log(abs(Ψ(x))) # Ψ=exp(S)    
    #
    Ψ_in = Ψ(R) # needed for fermion sign checks
    ∇S  = 0.5*drift(R)  # drift = 2∇S
    dim = size(R,1) 
    N = size(R,2)
    h = 1e-6
    Rloc = copy(R)
    
    for i in 1:N
        for k in 1:dim
            Rloc[k,i] += h            
            # check it's allowed, fermion Ψ may change sign.
            if Ψ(Rloc)*Ψ_in<0
                Rloc = copy(R)
                println("sign change, skip")
                break
            end
            S_p = S(Rloc)
            Rloc[k,i] -= 2h
            if Ψ(Rloc)*Ψ_in<0
                Rloc = copy(R)
                println("sign change, skip")
                break 
            end
            S_m = S(Rloc)
            Rloc[k,i] += h # move back
            ∇S_num = (S_p-S_m)/(2h)
            diff = (∇S[k,i]-∇S_num)/∇S[k,i]
            #@printf("∇S check: analytical: %20.12f numerical: %20.12f  rel difference: %20.12e\n",∇S[k,i],∇S_num,diff)
            if abs(diff)>1e-4 || isnan(diff)
                println(k," ",i)
                @printf("∇S check: analytical: %20.12f numerical: %20.12f  rel difference: %20.12e\n",∇S[k,i],∇S_num,diff)
                error("∇S check failed")
            end
        end
    end
end


function init_stat(datasize ::Int64, blocksize ::Int64; numblocks::Int64=100)
    stat = t_Stat()
    stat.blocksize = blocksize
    stat.nblocks = 0    
    stat.sample = t_StatData(0, zeros(datasize), zeros(datasize), 0.0)
    # data blocks    
    stat.datablock = Vector{t_StatData}(undef, numblocks) # was []
    stat.finished = false
    return stat
end


function add_sample!(stat ::t_Stat, dat)
    stat.finished  = false
    @. stat.sample.data  += dat
    @. stat.sample.data2 += dat^2
    stat.sample.n += 1
    if stat.sample.n == stat.blocksize
        # one full block collected, move average to block data
        # input data σ^2:
        N = stat.sample.n
        d2 = sum(stat.sample.data2)/N
        d  = sum(stat.sample.data)/N
        input_σ2 = d2-d^2
        #
        stat.nblocks += 1
        if stat.nblocks > length(stat.datablock)
            old_size =  length(stat.datablock)
            resize!(stat.datablock, old_size + 100) # unintialized elements in the end
        end 
        stat.datablock[stat.nblocks] = t_StatData(0, stat.sample.data./N, zeros(length(stat.sample.data)), input_σ2)
        # old, see "was" in init_stat
        #push!(stat.datablock,t_StatData(0, stat.sample.data./N, zeros(length(stat.sample.data)), input_σ2))        
        @. stat.sample.data = 0
        @. stat.sample.data2 = 0
        stat.sample.n = 0
        stat.finished = true
    end
end


function get_stats(stat ::t_Stat)
    N = stat.nblocks
    if N==0
        println("get_stats: no data")
        return 0, 0, 0, 0
    end
    if length(stat.datablock[1].data)==1
        ave_1, std_1, input_σ2_1, N_1 = get_stats_1(stat ::t_Stat)
        return ave_1, std_1, input_σ2_1, N_1,  stat.datablock[N].data[1] 
    end
    ave = similar(stat.datablock[1].data) 
    ave2 = similar(ave)
    ave .= 0
    ave2 .= 0
    input_σ2 = 0.0
    for i = 1:N
        d = stat.datablock[i].data
        @. ave += d
        @. ave2 += d^2
        input_σ2 += stat.datablock[i].input_σ2 
    end    
    @. ave /= N
    @. ave2 /= N
    input_σ2 /= N
    var = copy(ave)
    var2 = copy(ave)
    std = copy(ave)
    @. var2 = abs(ave2 - ave^2)
    @. var = sqrt(var2)
    @. std = var/sqrt(N)    
    return ave, std, input_σ2, N, stat.datablock[N].data 
end

function get_stats_1(stat ::t_Stat)
    N = stat.nblocks
    if N==0
        println("get_stats: no data")
        return nothing
    end
    ave = 0.0 
    ave2 = 0.0
    input_σ2 = 0.0 
    for i = 1:N
        d = stat.datablock[i].data[1]
        ave += d
        ave2 += d^2
        input_σ2 += stat.datablock[i].input_σ2 
    end 
    ave /= N
    ave2 /= N
    input_σ2 /= N
    var2 = abs(ave2 - ave^2)
    var = sqrt(var2) 
    std = var/sqrt(N)
    return ave, std, input_σ2, N
end

function output_MCresult(value, error)
    if isapprox(error,0.0)
        acc =  floor(Int64, log10(1/1e-10)+2)        
    else
        acc = floor(Int64, log10(1/error)+2)
    end
    fmt = Printf.Format("%."*"$(acc)f"*" +/- "*"%."*"$(acc)f \n")
    Printf.format(stdout, fmt, value, error)
    # convert +/-  to 12345(8) notation
    err = trunc(Int64,round(error*10^acc))
    fmt = Printf.Format("%."*"$(acc)f"*"("*"$(err)"*")\n\n")
    Printf.format(stdout, fmt, value)
end

function main()

    Random.seed!(1234)

    # Main program 

    # initialize
    println("He atom VMC")
    println("===========")
    walker, vmc_params = init()

    # choose trial wave function parameters defined in Model_Heatom.jl
    trial = :energy_optimized_parameters
    par = get_wave_function_params(trial)
    wf_params = (par.α, par.α12, par.β) # tuple
    @show wf_params 

    # closures with fixed wf_params
    Ψ_par(x) = Ψ(x, wf_params)
    EL_par(x) = EL(x, wf_params)
    drift_par(x) = drift(x, wf_params)

    rrs = zeros(Float64, 3, 2)
    d_buf = zeros(Float64, 3)

    # thermalization
    println("thermalizing") 
    for i in 1:Ntherm
        vmc_step!(walker.R, vmc_params, Ψ_par, rrs, d_buf)
        adjust_step!(vmc_params)
    end    
    println("thermalization done")
    println("-"^20," extra checks ","-"^20)
    println("Checking numerically local energy EL against trial wave function Ψ to spot errors in derivatives")
    
    # Checks could be done more accurately using AD
    for i in 1:5
        vmc_step!(walker.R, vmc_params, Ψ_par, rrs, d_buf)
        walker.E = EL_par(walker.R)
        num_check_EL(walker.R, walker.E, Ψ_par, eval_V)
    end
    num_check_∇S(walker.R, Ψ_par, drift_par)
    println("checks passed")
    println("-"^55)
    println()
    println()

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
        vmc_step!(walker.R, vmc_params, Ψ_par, rrs, d_buf)
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
            @printf("VMC %15d E = %.10f <E> = %.10f ± %.10f\n",ivmc, walker.E,  E_ave, E_std)

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
