using Printf
using Distributions
using StaticArrays
using Statistics


# local modules:
push!(LOAD_PATH,".")
using Utilities
using QMC_Statistics


const dim = 3 :: Int64
const D = (1/2) :: Float64
const N = 2   :: Int64 # number of particles

const Eexact = (9/2*π^2) ::Float64  # 44.41321980
const ET = 50.0 :: Float64 # a guess
const κ = 20.0 :: Float64 # may need adjustment
  
const Nwx = 20000 :: Int64
const Nw_target = 5000 :: Int64


if length(ARGS)==1
    τ_in = parse(Float64, ARGS[1])
else
    τ_in = 0.001 # default τ
end
const τ = τ_in :: Float64
const blocksize = floor(Int64,1/(10*τ)) :: Int64 # just an approx
const accuracy_goal = 0.05

function phi_1(x ::MVector)
    return sin(π*x[1])* sin(π*x[2])* sin(π*x[3])
end

function phi_2(x ::MVector)
    return sin(2*π*x[1])* sin(π*x[2])* sin(π*x[3])
end
    
function phi_3(x ::MVector)
    return sin(3*π*x[1])* sin(π*x[2])* sin(π*x[3])
end

# take *only* the nodes from the wave function
function psi_plus(x ::MMatrix)
    r1 = x[:,1]
    r2 = x[:,2]
    # exact GS
    psi = phi_1(r1)*phi_2(r2)-phi_1(r2)*phi_2(r1)
    # not GS
    #psi = phi_1(r1)*phi_3(r2)-phi_1(r2)*phi_3(r1)
    return psi>0
end        

mutable struct Walker
    R     ::MMatrix{dim,N,Float64}
    alive ::Int64
end

# initialization 
function init() ::Array{Walker,1}
    println("init")    
    # 
    # Initialize walkers
    #
    Rinit = @MMatrix zeros(Float64,dim,N)
    walker = [Walker(copy(Rinit),0) for i in 1:Nwx] # R, alive=0
    println("generating $Nw_target walkers")
    # set Nw_target walkers
    for iw in 1:Nw_target
        ww = walker[iw]
        ww.alive = 1
        while true
            R = @MMatrix rand(dim,N)        # coordinates
            R .* 0.8
            R .+ 0.1
            if psi_plus(R)
                ww.R .= R
                break
            end
        end        
        
    end
        
    println("init done")
    walker
end

# Diffusion
function diffusion!(x ::MMatrix)
    η = reshape(rand(Normal(0,1),dim*N),(dim,N))
    @. x += sqrt(2*D*τ)*η
end

# Branching
function branch!(walker ::Vector{Walker}, copies ::MVector)
    alives = [iw for iw in 1:Nwx if walker[iw].alive==1]
    deads = [iw for iw in 1:Nwx if walker[iw].alive==0]
    idead = 1
    for iw in alives
        #if walker[iw].alive==0; continue; end  # not alive, skip
        if copies[iw]==1; continue; end        # one copy is already there
        if copies[iw]==0
            # walker dies
            walker[iw].alive = 0
            continue
        end
        # copies[iw]>1
        # copy the walker to empty slots 
        for inew in 1:copies[iw]-1   # with copies=3 this is 1,2
            iw2 = deads[idead]
            walker[iw2].alive = 1
            walker[iw2].R = copy(walker[iw].R)
            idead += 1
            if idead>length(deads)
                error("cannot branch, outaspace")
            end
        end
    end
    Nw = 0
    @inbounds for iw in 1:Nwx
        if walker[iw].alive==1
            Nw += 1
        end
    end
    if Nw>=Nwx
        println("Error No free walker slots available; too large τ or too small κ?")
        exit(1)
    end 
    return Nw
end


    
    

function main()
    #
    # Main program 
    # 
    # initialize

    walker = init()

    #filename = "E_2fermions_in_a_box_no_importance_sampling_julia_13.dat"
    # filename = "E_2fermions_in_a_box_no_importance_sampling_julia.dat"
    filename = "E.scrap"
    println("output: ",filename)
    rm(filename, force=true)
    rm("E.tmp", force=true)
    
    #
    # DMC
    #
    Nw = Nw_target
    ET = 50.0  # trial energy, to be updated
    println("ET = $ET")
    
    Ntherm = floor(Int64,1.0/τ) # start DMC measurement after excitations have died out; system specific !
    idmc = -Ntherm

    E_ave = ET
    nE = 1
    

    copies = @MVector zeros(Int64,Nwx)
      
    # init E measurement
    Estat = init_stat(1, blocksize)

    while true
        #
        # take a DMC step in each walker
        #
        alives = Set([iw for iw in 1:Nwx if walker[iw].alive==1])
        @inbounds for iw in 1:Nwx
            if walker[iw].alive==0; continue; end  # not alive, skip
            #
            x = walker[iw].R
            
            diffusion!(x)            
            
            if any(i->i<0,x) || any(i->i>1,x)
                # outside the box
                copies[iw] = 0
            else
                if psi_plus(x)
                    copies[iw] = floor(Int64, exp(τ*ET) + rand())    
                else
                    # not a positive walker
                    copies[iw] = 0
                end
            end
            
            
        end
        
        # Branching
        Nw = branch!(walker, copies)
       
        idmc +=1
        
        # update trial energy; on average Nw_target walkers
        # A too large factor κ will cause a bad feedback,
        # a too small will let number of walkers get too large
        ET = E_ave/nE + κ*log(Nw_target/Nw)

        E_ave += ET 
        nE += 1
        
        # add new energy data
        add_sample!(Estat, ET)

        if idmc<=0
            @printf("DMC %10d  τ = %.5e  <E> = %.10f ET = %.10f <E>_exact = %.10f  %6d Walkers \n",
                    idmc, τ, E_ave/nE, ET, Eexact, Nw)
        end
        
        open("E.tmp","a") do f
            println(f,E_ave/nE)
        end
        # block data
        # ==========
        # screen and file output
        if Estat.finished            
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)
            @printf("DMC %10d  τ = %.5e <E> = %.10f +/- %.10f ET = %.10f <E>_exact = %.10f  %6d Walkers \n",
                    idmc, τ, Eb_ave, Eb_std, ET, Eexact, Nw)
            
            if Nb>10 && Eb_std < accuracy_goal
                println("reached accuracy goal")
                println("result: τ  <E>  error")
                @printf("%g ",τ)
                open(filename,"a") do f
                    println(f,τ," ",Eb_ave," ",Eb_std," ",Eexact)
                end
                output_MCresult(Eb_ave, Eb_std)
                break
            end

        end
        
        # 
        if idmc==0
            println("THERMALIZATION ENDS")
            Estat = init_stat(1, blocksize)
            E_ave = ET
            nE = 1
        end

      
    end
        
end
    
@time main()
