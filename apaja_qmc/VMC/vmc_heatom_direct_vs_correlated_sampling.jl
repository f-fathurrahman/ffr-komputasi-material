#=
JULIA

Correlated sampling energy optimization

VMC on a helium atom with trial wave function φ_T
defined in module Model_Heatom.jl

The Hamiltonian is
H = -1/2* (nabla_1^2 + nabla_2^2) - 2/r1 - 2/r2 + 1/r12  
=#


#using Distributions
using Printf
using StaticArrays
using Statistics
using LaTeXStrings
using BenchmarkTools
using InteractiveUtils
using Plots

using Optimization, OptimizationNLopt, NLopt
using Random

# local modules:
push!(LOAD_PATH,".")
using Common: VMC_Params, WfParams
using VMCstep
using Model_Heatom
using Utilities
#
# Parameters
#
const Ntherm = 100 # thermalization steps
const Nw = 100000 # number of walkers used in correlated sampling

# preallocated arrays
Esamples = Vector{Float64}(undef, Nw)
ratios = Vector{Float64}(undef, Nw)

mutable struct Walker
    R     ::MMatrix{dim,N,Float64}
    ψ     ::Float64
    E     ::Float64
end

#Random.seed!(123414) # For testing only 

# globals: 
#vmc_params = VMC_Params(0,0,0.0)
#walker = Vector{Walker}(undef, 1) # dummy init

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


function generate_walkers(wf_params::Tuple{Float64, Float64, Float64}, vmc_params ::VMC_Params, walker ::Vector{Walker})
    #println("generating walkers with wf_params ",wf_params)
    # thermalization
    Ψ_par(x) = Ψ(x, wf_params) # closure
    EL_par(x) = EL(x, wf_paarms)
    @inbounds for i in 1:Ntherm
        vmc_step!(walker[1].R, vmc_params, Ψ_par)      
    end
    # generate Nw walkers
    @inbounds begin
        for iw = 1:Nw
            for i = 1:10 # run walker 1 for a while
                vmc_step!(walker[1].R, vmc_params, Ψ_par)
            end
            ψR = vmc_step!(walker[1].R, vmc_params, Ψ_par)
            walker[iw].R = copy(walker[1].R)
            walker[iw].ψ = ψR
            walker[iw].E = EL(walker[1].R, wf_params)
        end
    end
end

function get_Eσ(wf_params::Tuple{Float64, Float64, Float64}, correlated::Bool, vmc_params ::VMC_Params, walker ::Vector{Walker}
                ; Esamples = Esamples, ratios = ratios)
    Esamples .= 0.0 
    if correlated        
        ratios .= 0.0 
        @inbounds for iw = 1:Nw
            wratio = Ψ(walker[iw].R, wf_params)^2/walker[iw].ψ^2         
            Esamples[iw] = EL(walker[iw].R, wf_params)* wratio
            ratios[iw] = wratio
        end
        r_ave = mean(ratios)
        E_ave = mean(Esamples)/r_ave
    else
        generate_walkers(wf_params, vmc_params, walker)  #  always generate new walkers 
        @inbounds for iw = 1:Nw
            Esamples[iw] = EL(walker[iw].R, wf_params)
        end
        E_ave = mean(Esamples)
    end    
    σ = std(Esamples)/sqrt(Nw)
    println("wf parameters:")
    for par in wf_params
        @printf("%20.15f  ", par)
    end
    println()
    

    if correlated
        @printf("CORRELATED SAMPLING <E> = %20.15f ± %20.15f  <weight ratio> = %20.15e\n", E_ave, σ, r_ave)
    else
         @printf("ORDINARY SAMPLING <E> = %20.15f ± %20.15f\n", E_ave, σ)
    end
    E_ave, σ
end


function main()
    #
    # Main program 
    #    
    walker, vmc_params = init()
   
    println("CORRELATED SAMPLING")
    par = get_wave_function_params(:initial_parameters_for_optimization)
    # seach limits
    β_min = 0.0
    β_max = 0.5
    Nβs = 50
    βs = LinRange(β_min, β_max, Nβs)

    # Define optimizer
    n_parameters = 1
    opt =  Opt(:LN_BOBYQA, n_parameters) 
    lower_bounds!(opt, β_min)
    upper_bounds!(opt, β_max)
    xtol_rel!(opt, 1e-5)    # Relative tolerance on optimization parameters
    maxeval!(opt, 10000)

    # help arrays
    Es = Array{Float64, 1}(undef, Nβs)
    σs = Array{Float64, 1}(undef, Nβs)
    
    # Correlated sampling
    # -------------------
    println("CORRELATED SAMPLING")
    generate_walkers((par.α, par.α12, par.β), vmc_params, walker)  # generate walkers only *once*
    for (i,β) in enumerate(βs)
        Es[i], σs[i] = get_Eσ((par.α, par.α12, β), true, vmc_params, walker) # true for correlated sampling        
    end
    p = plot(βs, Es, yerror = σs, xlabel = L"β",ylabel = L"E(β)",ylimits=(-2.90,-2.84), framestyle=:box,
             label="Correlated sampling")
    plot!([par.β], seriestype="vline",linecolor="red", label="Walkers sampled with this parameter")
    display(p)

    println("\n\nOptimization of E(β) using algorithm bobyqa\n")

    # E_opt(β) closure
    # x is an array, x[1] is β, the second [1] picks E from tuple output E, σ
    E_opt = x -> get_Eσ((par.α, par.α12, x[1]), true, vmc_params, walker)[1] 
    
    println("start optimization from β = ",par.β, " E = ", E_opt(par.β))
    
    # Define objective 
    min_objective!(opt, (x, _) -> E_opt(x))

    (E_corr_min, wf_params_corr_opt, ret) = optimize(opt, [par.β])

    # Output results
    println("minimum at β: ", wf_params_corr_opt)
    println("Minimum value: ", E_corr_min)
    println("Return code: ", ret) # FORCED_STOP emans failed, probably error in function evalution
    
    plot!([wf_params_corr_opt[1]], seriestype="vline", linecolor="green", label="Minimum found by algorithm")
    display(p)
    println("press enter")
    readline() # just to make the plots appears long enough to see them
    file = "correlated_sampling_He-atom.pdf"
    println("output to file $file")
    savefig(p, file)

    # Ordinary sampling
    # -----------------
    println("ORDINARY SAMPLING")
    
    correlated = false
    for (i,β) in enumerate(βs)
        Es[i], σs[i] = get_Eσ((par.α, par.α12, β), correlated, vmc_params, walker)        
    end
    p = plot(βs, Es, yerror = σs, xlabel = L"β",ylabel = L"E(β)", ylimits=(-2.90,-2.84), framestyle=:box,
             label="Ordinary sampling")
       
    println("\n\nOptimization of E(β) using algorithm bobyqa\n")
    
    E_opt = x -> get_Eσ((par.α, par.α12, x[1]), false, vmc_params, walker)[1]  # false for non-correlated samples 
    min_objective!(opt, (x, _) -> E_opt(x))
    println("start optimization from β = ",par.β, " E = ", E_opt(par.β))
    (E_min, wf_params_opt, ret) = optimize(opt, [par.β])
    println("minimum at β: ", wf_params_opt)
    println("Minimum value: ", E_min)
    println("Return code: ", ret)
    println("Notice, that the optimization algorithm evaluates E(β) wherever it wants,")
    println("and these (β, E(β)) are not those shows in the plot.") 
    plot!([par.β], seriestype="vline",linecolor="red", label="Walkers sampled with this parameter")
    plot!([wf_params_corr_opt[1]], seriestype="vline", linecolor="green", label="Minimum found by algorithm")
    display(p)
    println("press enter")
    readline() # just to make the plots appears long enough to see them
    file = "ordinary_sampling_He-atom.pdf"
    println("output to file $file")
    savefig(p, file)

    
end


@time main()
