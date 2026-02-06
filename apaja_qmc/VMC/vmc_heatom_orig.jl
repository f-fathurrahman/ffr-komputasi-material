#=
JULIA

VMC on a helium atom with trial wave function φ_T
defined in module Model_Heatom

The Hamiltonian is
H = -1/2* (nabla_1^2 + nabla_2^2) - 2/r1 - 2/r2 + 1/r12  
=#

import LinearAlgebra: norm
import Random
using Printf
using StaticArrays
using Statistics

# local modules:
push!(LOAD_PATH, ".")

using Common: VMC_Params
using VMCstep
using Model_Heatom
using Utilities
using QMC_Statistics
#
# Parameters
#
const blocksize = 1000000   # data block size
const Ntherm = 100         # thermalization steps
const accuracy_goal = 2e-4

mutable struct Walker
    R     ::MMatrix{dim,N,Float64}
    ψ     ::Float64
    E     ::Float64
end


# initialization 
function init()     
    # 
    # Initialize walker
    #
    R = @MMatrix rand(dim,N)    # coordinates
    walker = Walker(R, 0.0, 0.0) # R, ψ=0, E=0    
    vmc_params = VMC_Params(0,0,3.1)
    walker, vmc_params
end

       
function main()

    Random.seed!(1234)

    #
    # Main program 
    # 
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
    # thermalization
    println("thermalizing")    
    for i in 1:Ntherm
        vmc_step!(walker.R, vmc_params, Ψ_par)
        adjust_step!(vmc_params)
    end    
    println("thermalization done")

    println("-"^20," extra checks ","-"^20)
    println("Checking numerically local energy EL against trial wave function Ψ to spot errors in derivatives")
    # Checks could be done more accurately using AD
    for i in 1:5
        vmc_step!(walker.R, vmc_params, Ψ_par)
        walker.E = EL_par(walker.R)
        num_check_EL(walker.R, walker.E, Ψ_par, V)
    end
    num_check_∇S(walker.R, Ψ_par, drift_par)
    println("checks passed")
    println("-"^55)
    println()
    println()
    
    # saving to file is slow
    #filename = "E_heatom_VMC"
    #println("Storing energies to file $filename (deleting old file)")
    #rm(filename, force = true)
    #
    # VMC
    #
    # init E measurement (QMC_Statistics)
    Estat = init_stat(1, blocksize)
    #
    ivmc = 0
    while true # until accuracy_goal is reached
        vmc_step!(walker.R, vmc_params, Ψ_par)        
        walker.E = EL_par(walker.R)

        # add new energy data
        add_sample!(Estat, walker.E)

        ivmc +=1
        if ivmc%10 == 0                      
            adjust_step!(vmc_params)
        end
        # saving to file is slow
        #open(filename,"a") do f
        #    println(f,ivmc," ",E)
        #end
        #
        # output when a block is full
        if Estat.finished            
            E_ave, E_std, E_inputvar2, Nb = get_stats(Estat)
            @printf("VMC %15d E = %.10f <E> = %.10f ± %.10f\n",ivmc, walker.E,  E_ave, E_std)

            if Nb>10 && E_std < accuracy_goal
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
end

# main()
