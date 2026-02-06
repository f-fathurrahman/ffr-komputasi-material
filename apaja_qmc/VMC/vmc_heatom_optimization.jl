#=
JULIA

Correlated sampling energy optimization

VMC on a helium atom with trial wave function φ_T
defined in module Model_Heatom.jl

The Hamiltonian is
H = -1/2* (nabla_1^2 + nabla_2^2) - 2/r1 - 2/r2 + 1/r12  
=#


using Printf
using StaticArrays
using Statistics
using LaTeXStrings

using Plots

using Optimization, OptimizationNLopt, NLopt
using Random
#using PRIMA # May 2025: requires executable stack, which is disabled in many systems for security reasons

# local modules:
push!(LOAD_PATH,".")

using Common
using VMCstep
using Model_Heatom
using Utilities

#
# Parameters
#
const Ntherm = 100 # thermalization steps
const Nw = 1_000_000 # number of walkers used in correlated sampling

mutable struct Walker
    R     ::MMatrix{dim,N,Float64}
    ψ     ::Float64
    E     ::Float64
end

#Random.seed!(123414) # For testing only

 
# initialization 
function init() 
    println("init ",Nw," walkers")    
    # 
    # Initialize walkers
    #
    R = @MMatrix rand(dim,N)    # coordinates
    walker = [Walker(R, 0.0, 0.0) for i in 1:Nw] # R, ψ=0, E=0    
    vmc_params = VMC_Params(0,0,3.1)
    walker, vmc_params
end


# Correlated sampling functions
# -----------------------------
function generate_walkers(wf_params::Vector{Float64}, vmc_params ::VMC_Params, walker ::Vector{Walker})   
    # thermalization
    Ψ_par(x) = Ψ(x, wf_params) # closure
    for i in 1:Ntherm
        vmc_step!(walker[1].R, vmc_params, Ψ_par)      
    end
    println("generating walkers")
    # generate Nw walkers
    Eave = 0.0
    for iw = 1:Nw
        for i = 1:10 # run walker 1 for a while
            vmc_step!(walker[1].R, vmc_params, Ψ_par)            
        end
        ψ = vmc_step!(walker[1].R, vmc_params, Ψ_par)            
        walker[iw].R = copy(walker[1].R)
        walker[iw].ψ = ψ
        walker[iw].E = EL(walker[1].R, wf_params)
        Eave += walker[iw].E
    end
    @printf("initial <E> = %15.5f\n", Eave/Nw) 
    println("done")
end
    
function get_correlated_Eσ(wf_params::Vector{Float64}, walker::Vector{Walker})
    Esamples = Vector{Float64}(undef, Nw)
    ratios = Vector{Float64}(undef, Nw)
    for iw = 1:Nw
        wratio = Ψ(walker[iw].R, wf_params)^2/walker[iw].ψ^2         
        Esamples[iw] = EL(walker[iw].R, wf_params)* wratio
        ratios[iw] = wratio
    end
    r_ave = mean(ratios)
    E_ave = mean(Esamples)/r_ave
    
    σ = std(Esamples)/sqrt(Nw)
    for par in wf_params
        @printf("wf parameter %20.15f\n",par)
    end
    @printf("correlated energy = %20.15f ± %20.15f  <weight ratio> = %20.15e\n", E_ave, σ, r_ave)
    E_ave, σ
end

# driver for bobyqa
function make_correlated_E_driver(walker)
    wf_params -> get_correlated_Eσ(wf_params, walker)[1] # use only 1st of the tuple
end

function get_correlated_σ2(wf_params, walker::Vector{Walker})
    Eguess = 0.0
    for iw = 1:Nw
        Eguess += walker[iw].E
    end
    Eguess /= Nw
    
    σ2samples = Array{Float64, 1}(undef, Nw)
    ratios = Array{Float64, 1}(undef, Nw)    
    for iw = 1:Nw
        wratio = Ψ(walker[iw].R, wf_params)^2/walker[iw].ψ^2         
        σ2samples[iw] = (EL(walker[iw].R, wf_params)-Eguess)^2 * wratio
        ratios[iw] = wratio
    end
    r_ave = mean(ratios)
    σ2_ave = mean(σ2samples)/r_ave
    σ = std(σ2samples)/sqrt(Nw)
    for par in wf_params
        @printf("wf parameter %20.15f\n",par)
    end
    @printf("σ^2 = %20.15f ± %20.15f  <φ_T ratio> = %20.15e\n", σ2_ave, σ , r_ave)
    σ2_ave, σ
end

# driver for bobyqa
function make_correlated_σ_driver(walker)
    wf_params -> get_correlated_σ2(wf_params, walker)[1] # use only 1st of the tuple
end



function main()
    #
    # Main program 
    #
   
    walker, vmc_params =  init()
   
    println("CORRELATED SAMPLING")
    par = get_wave_function_params(:initial_parameters_for_optimization)
    α = par.α
    β = par.β
    α12 = par.α12
    println("\nOptimization of E\n")
    # start with values from Model
    u0 = [α, α12, β]   
    wf_params = copy(u0)  

    println("start with parameters = ", wf_params)
    generate_walkers(wf_params, vmc_params, walker)  # generate walkers only *once*

    E_u0, σ_u0 = get_correlated_Eσ(wf_params, walker)

    
    @printf("initial correlated energy %15.5f variance %15.5f\n", E_u0, σ_u0)
    # limits
    # optimize all parameters
    xl = [0.0, 0.0, 0.0]
    xu = [20.0, 20.0, 20.0]
    
    # optimize only α and β:
    #xl = [0.0, α12, 0.0]
    #xu = [3.0, α12+1e-5, 1.0]
    # optimize only β:
    #xl = [α, α12, 0.0]
    #xu = [α+1e-5, α12 + 1e-5, 1.0]  # algorithm needs finite range


    # Define optimizer
    n_parameters = 3
    opt =  Opt(:LN_BOBYQA, n_parameters) 
    lower_bounds!(opt, xl)
    upper_bounds!(opt, xu)
    xtol_rel!(opt, 1e-5)    # Relative tolerance on optimization parameters
    maxeval!(opt, 100000)   # Maximum number of evaluations
    
    # Define objective
    E_opt_correlated =  make_correlated_E_driver(walker)    
    min_objective!(opt, (x, _) -> E_opt_correlated(x))

    (E_ene_min, wf_params_ene_min, ret) = optimize(opt, u0)

    # Output results
    println("Optimized parameters: ", wf_params_ene_min)
    println("Minimum value: ", E_ene_min)
    println("Return code: ", ret)
        
    # Variance optimization
    # ---------------------
    # use same walkers as in energy optimization
       

    println("\n\nOptimization of σ^2\n")
    
    u0 = [α, α12, β] # start with values from Model
    wf_params = copy(u0)
    println("start with parematers = ", wf_params)

    # Define objective
    σ2_opt_correlated =  make_correlated_σ_driver(walker)    
    min_objective!(opt, (x, _) -> σ2_opt_correlated(x))

    (σ2_min, wf_params_var_min, ret) = optimize(opt, u0)

    # Output results
    println("Optimized parameters: ", wf_params_var_min)
    println("Minimum value: ", σ2_min)
    println("Return code: ", ret)
    
    println("\n\nRESULTS:")
    println("Energy minimum:", E_ene_min," " , wf_params_ene_min )
    println("Variance minimum:", σ2_min, " ", wf_params_var_min)
    
end


@time main()
