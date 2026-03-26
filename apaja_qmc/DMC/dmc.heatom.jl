#=
JULIA

 DMC on a helium atom with trial wave function 
 φ_T = e^{-α*r1}e^{-α*r2} ; r1,r2 distance of els from proton
 one el in spin-up and the other in spin-down state
  
Measures only total energy
Hamiltonian in a.u. (Z=2) :

H = -1/2* (nabla_1^2 + nabla_2^2) - 2/r1 - 2/r2 + 1/r12

run with τ=0.001 (default) using
  julia dmc.heatom
run with τ = 0.01 using 
  julia dmc.heatom 0.01 
=#


import LinearAlgebra: norm

using Distributions
#using LoopVectorization
using Printf
using StaticArrays
using Statistics
using BenchmarkTools


# local modules:
push!(LOAD_PATH,".")
using VMCstep
using Utilities: dist1, dist2, num_check_EL, num_check_∇S, metro, output_MCresult
using QMC_Statistics
using Common: VMC_Params
#
# Parameters
#
const order = 1   # 0, 1 or 2
println("He atom DMC using algorithm of order ",order)

# He atom parameters
const D = 0.5                # hbar^2/(2m) in a.u.
const N = 2                  # number of electrons
const Z = 2                  # change number  
const dim = 3                # dimension           
const Eexact = -2.903724377  # accurate ground-state energy


# VMC and DMC parameters

# imaginary time step, either default or command line argument 
t = 0.001         #  default
# read τ from a command line argument
if length(ARGS)!=0
    t = parse(Float64, ARGS[1]) # convert string to float
end
const τ::Float64 = t 

const blocksize = floor(Int64,100/τ) # DMC data block size

const NVMC = 1000        # number of VMC steps
const NDMC = 10000000    # number of DMC steps
const Nwx = 20000        # max number of walkers
const Nw_target = 1000   # on average this many walkers
const κ = 0.1            # how zealously DMC tries to keep Nw_target walkers  

# ==================================
# Trial wf
# ========
# parameters
const α = 2.0 # 1.847529
const β = 0.159321
const γ = 1.0
const α12 = 0.5 # 0.359070

@inline function get_unitvecs(R::MMatrix{dim, N, Float64}) where {dim, N}
    vecr1 = R[:, 1]
    vecr2 = R[:, 2]
    vecr12 =  vecr1 - vecr2
    r1 = dist1(R, 1)
    r2 = dist1(R, 2)
    r12 = dist2(R, 1, 2)
    vecr1/r1, vecr2/r2, vecr12/r12, r1, r2, r12
end

function V(R::MMatrix{dim, N, Float64}) where {dim, N}
    r1 = dist1(R, 1)
    r2 = dist1(R, 2)
    r12 = dist2(R, 1, 2) 
    - Z/r1 - Z/r2 + 1/r12
end

@inline function φ_T(R::MMatrix{dim, N, Float64}) where {dim, N}
    r1 = dist1(R, 1)
    r2 = dist1(R, 2)
    r12 = dist2(R, 1, 2)
    b = 1+β*r12
    exp(-α*(r1+r2) + α12*r12/b)
end

const ∇S_buf = MMatrix{3, 2, Float64}(undef)
# F=2∇S
@inline function F(R::MMatrix{dim, N, Float64}) where {dim, N}
    hatr1, hatr2, hatr12, r1, r2, r12 =  get_unitvecs(R)
    b = 1+β*r12    
    term12 = α12/b^2*hatr12
    ∇S = ∇S_buf 
    @inbounds begin        
        ∇S[:, 1] = -α * hatr1 + term12
        ∇S[:, 2] = -α * hatr2 - term12
    end
    2∇S
end

@inline function EL(R::MMatrix{dim, N, Float64}) where {dim, N}
    hatr1, hatr2, hatr12, r1, r2, r12 = get_unitvecs(R)
    b = 1+β*r12
    ∇S = 0.5*F(R)
    ∇2S = -α*2/r1 -α*2/r2  + 4*α12/b^3*1/r12 
    EL = -0.5*(sum(∇S.^2) + ∇2S) - Z/r1 - Z/r2 + 1/r12
end

# ==================================


# local Walker struct
mutable struct Walker
    R ::MMatrix{dim, N, Float64}
    E ::Float64
end

# initialization 
function init()
    println("init")    
    # 
    # Initialize walkers
    #
    walker = [Walker(MMatrix{dim, N, Float64}(undef), 0.0) for i in 1:Nwx]
    println("generating $Nw_target walkers")
    alives = falses(Nwx)  # BitVector
    # set Nw_target walkers
    for iw in 1:Nw_target
        R = @MMatrix rand(dim, N)   # coordinates        
        walker[iw].R .= R
        walker[iw].E = EL(R)
        alives[iw] = true
    end
    vmc_params = VMC_Params(0, 0, 3.1) # initial VMC step 3.1    
    println("init done")
    walker, alives, vmc_params
end

# branch for the local Walker struct
function branch!(walker ::Vector{Walker}, copies ::MVector, alives ::BitVector)
    # list initially dead walkers that can we can branch to
    deads = .!alives
    alive_walkers = findall(alives)
    @inbounds for iw in alive_walkers
        if copies[iw]==1
            continue        # one copy is already there
        end
        if copies[iw]==0
            # walker dies
            alives[iw] = false
            deads[iw] = true
            continue
        end
        # copies[iw]>1
        # copy the walker to empty slots 
        @inbounds for inew in 1:copies[iw]-1  
            # use free slot
            iw2 = findfirst(deads)
            if iw2 == nothing
                error("No free walkers available; too large τ ?")
            end
            walker[iw2].R .= walker[iw].R
            walker[iw2].E =  walker[iw].E            
            alives[iw2] = true
            deads[iw2] = false
        end
    end
    return nothing
end

# one diffusion+drift step in DMCstep = vmc_params.step
    
@inline function diffusion_drift_step!(x ::MMatrix{dim, N, Float64},
                                       y ::MMatrix{dim, N, Float64},
                                       z ::MMatrix{dim, N, Float64}) where {dim, N}   
    # help arrays y and z
    if order==1
        # diffusion(τ)+drift(τ)
        η = randn(dim, N)
        Fx = F(x) 
        @. x += sqrt(2*D*τ)*η + D*τ*Fx 
           
    elseif order==2
        #
        # diffusion(τ/2)+drift(τ)+diffusion(τ/2)
        #
        # step 1)
        η = randn(dim, N)
        @. y = x + sqrt(D*τ)*η    
        Fy = F(y)
        # step 2)
        @. z = y + D*τ/2*Fy        
        Fz = F(z)
        # step 3)
        η = randn(dim, N)
        @. x = y + D*τ*Fz + sqrt(D*τ)*η
    elseif order==3
        dt = 2*D*τ
        eta = randn(dim, N)
        @. x = x + sqrt(dt/2)*eta
        F1 = F(x)
        @. y = x + dt/2*F1/2 
        F2 = F(y)
        eta = randn(dim, N)
        @. x += dt*F2/2 + sqrt(dt/2)*eta
    end

    
end



function lnG(xp ::MMatrix{dim, N, Float64}, x ::MMatrix{dim, N, Float64}) where {dim, N}        
    # only parts not symmetric in x'<->x
    # G(x'<-x,τ) = exp(-(x'-x-dτF(x))^2/(4Dτ)
    # 
    if order != 1
        error("lnG only for 1st order code")
    end
    lnG = -norm(xp - x - D*τ*F(x))^2 /(4*D*τ) 
end

        
function main()
    #
    # Main program 
    # 
    # initialize
    walker, alives, vmc_params = init() 
    # thermalization
    println("thermalizing")
    Nw = Nw_target # for VMC
    @inbounds for i in 1:10
        @inbounds for iw in 1:Nw
            vmc_step!(walker[iw].R, vmc_params, φ_T)
        end
    end

    println("thermalization done")

    println("Checking numerically EL and drift against φ_T")
    for i in 1:5
        vmc_step!(walker[1].R, vmc_params, φ_T)
        R = walker[1].R
        E = EL(R)
        num_check_EL(R, E, φ_T, V)
    end

    num_check_∇S(walker[1].R, φ_T, F)
    println("checks passed")

    filename = string("E_heatom_",order,"_tau=",τ)
    #filename = string("E_heatom_noacceptreject",order,"_tau=",τ)
    println("output: ",filename)
    # init file output
    open(filename,"w") do f
        println(f," ")
    end
    #
    # VMC
    #
    E_ave = 0
    nE  = 0
    @inbounds for ivmc in 1:NVMC
        E = 0.0 
        @inbounds for iw in 1:Nw         
            vmc_step!(walker[iw].R, vmc_params, φ_T)
            walker[iw].E = EL(walker[iw].R)
            E += walker[iw].E 
        end
        E_ave += E/Nw        
        nE +=1
        
        if ivmc%10 == 0            
            @printf("VMC  E = %.10f <E> = %.10f\n",E/Nw,E_ave/nE)            
        end
        adjust_step!(vmc_params)
    end
    
    #
    # DMC
    #
    ET = E_ave/nE  # trial energy, to be updated
    println("ET = $ET")
    
    Ntherm = floor(Int64, 1.0/τ) # start DMC measurement after excitations have died out; system specific !
    idmc = -Ntherm

    E_ave = 0
    nE = 0
        
    Rold = @MMatrix zeros(Float64, dim, N)

    copies = @MVector zeros(Int64, Nwx)
    Ntry::Int64 = 0 
    Nacc::Int64 = 0 

    help1 = @MMatrix zeros(Float64, dim, N)
    help2 = @MMatrix zeros(Float64, dim, N)
   
    # init E measurement
    Estat = init_stat(1, blocksize)

    while true
        #
        # take a DMC step in each walker
        #
        
        @inbounds for iw in findall(alives)
            ELold::Float64 = walker[iw].E   
            copy!(Rold, walker[iw].R) 
            Rnew = walker[iw].R # alias
            #
            diffusion_drift_step!(Rnew, help1, help2)
            # 
            
            ELnew::Float64 = EL(walker[iw].R)    
            
            accept = true
            
            if order==1
                # metropolis-Hastings                
                Wnew::Float64 = lnG(Rold, Rnew) + 2*log(φ_T(Rnew))  # lnG(new->old) + lnpsi2(new)
                Wold::Float64 = lnG(Rnew, Rold) + 2*log(φ_T(Rold))  # lnG(old->new) + lnpsi2(old)
                accept = metro(Wold,Wnew)
            end
            Ntry +=1
            if accept
                Nacc +=1
                walker[iw].E = ELnew
            else
                ELnew = ELold
                walker[iw].E = ELold
                Rnew .= Rold                
            end

            weight::Float64   = exp(-τ*(0.5*(ELold+ELnew)-ET))  # symmetric under R<->R'
            copies[iw] = floor(Int64, weight + rand()) 
        end
        # Branching
        branch!(walker, copies, alives)
        
        # collect energy from walkers
        # note: do this *after* drift, diffusion and branching ! 
        E = 0.0  
        @inbounds for iw in findall(alives)
            E  +=  walker[iw].E
        end

        
        Nw = sum(alives)             
        if Nw == Nwx
            println("Error: hit max number of walkers")
            exit(1)
        end
        E /= Nw
        E_ave += E 
        nE += 1
        # add new energy data
        add_sample!(Estat, E)
        idmc += 1
        
        # update trial energy; on average Nw_target walkers
        # A too large factor κ will cause a bad feedback,
        # a too small will let number of walkers get too large
        ET = E_ave/nE + κ*log(Nw_target/Nw)       

        if idmc<=0
            @printf("DMC %10d  E = %.10f <E> = %.10f ET = %.10f Eexact = %.10f  %6d Walkers \n",
                    idmc, E, E_ave/nE, ET, Eexact, Nw)
            if order==1;@printf("DMC acceptance = %.5f %%\n",Nacc*100.0/Ntry);end
        end

        # block data
        # ==========
        # screen and file output
        if Estat.finished
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)
            @printf("DMC %10d E %.10f <E> = %.10f +/- %.10f ET = %.10f Eexact = %.10f  %6d Walkers \n",
                    idmc, E ,Eb_ave, Eb_std, ET, Eexact, Nw)
            open(filename,"a") do f
                println(f,τ," ",Eb_ave," ",Eb_std," ",Eexact)
            end
        end
        
        # 
        if idmc==0
            println("THERMALIZATION ENDS")
            Estat = init_stat(1, blocksize)
            E_ave = 0
            nE = 0
        end

        if idmc == NDMC
            println("Trial wf parameters:")
            @show(α)
            @show(α12)
            @show(β)
            println("result: τ  <E>  error")
            @printf("%g ",τ)
            output_MCresult(Eb_ave, Eb_std)
            break
        end
    end
        
end
    
@time main()

