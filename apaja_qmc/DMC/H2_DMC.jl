
#
# Equilibrium of H2 molecule
# trial wave function almost similar to Traynor, Anderson and Boghosian (TAB) 
#
#
# Run, for example, as
# (i) optimization of parameters:
#  julia H2_DMC.jl optimize=1
# (ii) DMC, using time step 0.001 
#  julia H2_DMC.jl tau=0.001  
#
#  default Req=1.4 

#
#
# my results:
# for Params(1.282421875, 0.398046875, 1.3885, 11.24375)
# VMC:  <E> =  -1.14651 +/- 0.00093
#
using Printf
using StaticArrays
import Distributions: Normal
using Statistics
using ForwardDiff
using LinearAlgebra: diag, tr
using Random 
using Optim

@inline norm(x) = sqrt(sum(abs2,x))

# local modules
push!(LOAD_PATH,".")
using Common: VMC_Params, Walker, branch!
using LinearOptimization
using QMC_Statistics
using Utilities: argparse, num_check_EL, num_check_∇S, output_MCresult
using VMCstep


# parameters
const dim=3
const N=4
# hbar^2/(2m) for electrons and protons
const D = [1/2, 1/(2*1836.152673426)]

# Threads and RNG
println("Using $(Threads.nthreads()) threads")
# seed rng's using time in milliseconds
base_number = floor(Int, time() * 1000) 
base_number = 12345
rngs = [MersenneTwister(base_number + i) for i in 1:Threads.nthreads()]



# Wave function parameters
mutable struct Params{T<:Real}
        ζ :: T  # AO orbital parameter
        b :: T  # e-e correlation parameter
        c :: T  # p-p correlation parameter
        d :: T  # p-p correlation parameter, harmonic coeff  
end

function update_params(p ::Params, Δp ::MVector)    
    return Params(p.ζ+Δp[1], p.b+Δp[2], p.c+Δp[3], p.d+Δp[4])
end

# Command line parameters
# either 
#   optimize=1
# or
#   tau= (give imaginary time step) 
# 
possible_args=["optimize","tau"]
arg_dict = argparse(possible_args)
optimize = get(arg_dict, "optimize", 0)

if optimize==1
    const outfile_opt = "H2_TAB_opt.dat"
    rm(outfile_opt, force=true)    
else

    # use command line values or defaults
    τ_in = get(arg_dict, "tau", 0.0)
    if τ_in <1e-10
        error("DMC : set tau in command line")
    end
    
    # DMC parameters
    const τ = τ_in :: Float64
    const κ = 20.0 :: Float64 # may need adjustment
    const accuracy_goal = 0.00001
    const Nwx = 20000 :: Int64
    const Nw_target = 5000 :: Int64
    const blocksize = floor(Int64,1/τ) :: Int64 # just an approx
    const Ntherm = floor(Int64,10/τ)
    
    @show blocksize
    @show Ntherm
    
    const Efile = "E.H2_TAB_"*string(τ)
    rm(Efile, force=true)
    const Eallfile = "E.tmp_"*string(τ)
    rm(Eallfile, force=true)
  
end

# measurements
# ============

const measurement_list = [:Req, :Vpot, :dens_profile]
# sizes of measured data
const data_size = Dict(:Req=>1, :Vpot=>1, :dens_profile=>301)
# frequencies of measurements, expensive measurements infrequent
const measure_freq = Dict(:Req=>5, :Vpot=>5, :dens_profile=>10)

dens_grid = range(-3.0, 4.0, length=data_size[:dens_profile])

# Actual measurements
function measure_Req(walker ::Vector{Walker})
    ave = 0.0
    n = 0
    for w in walker
        if w.alive==1
            value = norm(w.R[:,3]-w.R[:,4])
            w.measurements[:Req] = value
            ave += value
            n += 1
        end
    end
    return ave/n
end

function measure_Vpot(walker ::Vector{Walker})
    ave = 0.0
    n = 0
    for w in walker
        if w.alive==1
            value = V(w.R)
            w.measurements[:Vpot] = value
            ave += value
            n += 1
        end
    end
    return ave/n
end

# electron density n(x) in proton-proton cylinder of radius rcyl
function measure_dens_profile(walker ::Vector{Walker})
    dens_profile = zeros(data_size[:dens_profile])
    n = 0
    rcyl = 0.1
    for w in walker
        if w.alive==1
            n += 1 
            r3 = w.R[:,3]
            r4 = w.R[:,4]
            r34 = r3-r4
            d34 = norm(r34)
            r34 = r34/d34
            # proton r3 is at 0.0 and proton r4 is at 1.0
            for i in 1:2
                r = (w.R[:,i] - r3) ./ d34 # set r3 at origin, scale distances with d34
                # project on cylinder axis
                r_proj = sum(r .* r34) .* r34
                # distance from cylinder axis
                dist = norm(r-r_proj)
                if dist<rcyl
                    # electron in the cylinder                    
                    # bin values
                    if sum(r_proj .* r34)<0.0
                        x = norm(r_proj)
                    else
                        x = -norm(r_proj)
                    end
                    bin = findfirst(b -> b > x, dens_grid)
                    if bin != nothing
                        dens_profile[bin] += 1
                    end
                end
            end
        end
    end
    bin_width = (dens_grid[end] - dens_grid[1]) / length(dens_grid)
    return dens_profile./ (n*bin_width)
end

function save_density_profile(prefix ::String, ave ::Vector{Float64}, std ::Vector{Float64})
    measure = :dens_profile
    filename = prefix*"."*string(measure)*".dat"
    println("saving density profile to $filename")
    open(filename,"w") do f
        for i in 1:data_size[measure]
            println(f, dens_grid[i]," ",ave[i]," " ,std[i])
        end
    end
end


measurement_functions = Dict(
    :Req  => measure_Req,
    :Vpot => measure_Vpot,
    :dens_profile => measure_dens_profile
)

# block data collections
function init_stat_collection!(measure ::Symbol, ndata ::Int64, blocksize ::Int64,  stat_collection ::Dict)
    stat_collection[measure] = init_stat(ndata, blocksize)
end

# all measurements
function init_stat_collection!(blocksize ::Int64)
    stat_collection = Dict{Symbol, Any}()
    for measure in measurement_list        
        init_stat_collection!(measure, data_size[measure], blocksize, stat_collection ::Dict)
    end
    return stat_collection
end

function do_measurements(walker ::Vector{Walker}, steps ::Int64, stat_collection ::Dict)
    for measure in measurement_list
        if steps%measure_freq[measure] == 0  # some quantities are measured less frequently
            measure_func = measurement_functions[measure]
            ave = measure_func(walker) # sets all w.measurements[measure] 
            add_sample!(stat_collection[measure], ave)
        end       
    end
end


# =================================



function φ_T(r ::MMatrix, p ::Params) 
    ζ, b, c, d = p.ζ, p.b, p.c, p.d
    r1, r2, r3, r4 = r[:,1], r[:,2], r[:,3], r[:,4]
    # e-p correlations
    r13  = norm(r1-r3)
    r14  = norm(r1-r4)
    r23  = norm(r2-r3)
    r24  = norm(r2-r4)
    wf = (exp(-ζ*r13)+exp(-ζ*r14))*(exp(-ζ*r23)+exp(-ζ*r24))
    # e-e correlations
    r12 = norm(r1-r2)
    wf *= exp(r12/(2*(1+b*r12))) # not quite what TAB used
    # p-p correlations
    r34 = norm(r3-r4)
    wf *= exp(-d*(r34-c)^2)
end

# local kinetic energy TL =-D  ∇^2 ϕ(r) / ϕ(r)
function TL_hessian_version(r ::MMatrix, ϕ ::Function)
    # for possible thread safety issues:
    #config = ForwardDiff.HessianConfig(ϕ, r)
    #hess =  ForwardDiff.hessian(ϕ, r, config) # size is (12,12)
    #
    hess = ForwardDiff.hessian(ϕ, r) # size is (12,12)
    dN = size(hess,1)
    dNhalf = Int(dN/2)
    # electrons and protons have different D:=hbar^2/(2m)   
    D∇2ϕ_r = D[1]*tr(hess[1:dNhalf, 1:dNhalf]) + D[2]*tr(hess[dNhalf+1:dN, dNhalf+1:dN])
    return -D∇2ϕ_r/ϕ(r) 
end

# non-hessian version
function TL(r ::MMatrix, ϕ ::Function)
    Dsum_∇2ϕ = 0.0
    dim, N = size(r)
    Nhalf = Int(N/2)
    for i in 1:N
        if i <= Nhalf
            Di = D[1]
        else
            Di = D[2]
        end
        for k in 1:dim
            function f(t)
                r_copy = similar(r, typeof(t)) 
                r_copy .= r  
                r_copy[k, i] = t                
                return ϕ(r_copy)      
            end
            ∇2  = ForwardDiff.derivative(t -> ForwardDiff.derivative(f, t), r[k, i])
            Dsum_∇2ϕ += Di*∇2
        end
    end
    return -Dsum_∇2ϕ/ϕ(r)
end


# potential energy
@inline function V(r ::MMatrix{dim,N,Float64}) where {dim, N}
    r1, r2, r3, r4 = r[:,1], r[:,2], r[:,3], r[:,4]
    # p-p
    r34  = norm(r3-r4)
    V = 1/r34
    # e-p 
    r13  = norm(r1-r3)
    r14  = norm(r1-r4)
    r23  = norm(r2-r3)
    r24  = norm(r2-r4)
    V += -1/r13 - 1/r14 - 1/r23 - 1/r24
    # e-e
    r12 = norm(r1 -r2) 
    V += 1/r12
    V
end


# local energy
function EL(r ::MMatrix, ϕ ::Function, params ::Params)
    return TL(r, x -> ϕ(x, params)) + V(r)    
end


function run_VMC!(r ::MMatrix{dim,N,Float64}, ϕ::Function,
                  vmc_params ::Vector{VMC_Params}, params ::Params;
                  tol = 1e-2, verbose = false )  where {dim, N}  

    @inbounds for i in 1:1000
        vmc_step!(r, vmc_params,  x -> ϕ(x, params))
        if i%10 == 0
            for p in vmc_params
                adjust_step!(p)
            end
        end
    end
    if verbose
        for p in vmc_params
            println("VMC       step = $(p.step)")
            println("VMC acceptance = $(p.Naccept*100.0/p.Ntry)")
        end
    end
    #
    # VMC
    #
    Estat = init_stat(1, 100)
    nn = 1000
    while true
        E_ave = 0.0 
        @inbounds for i in 1:nn
            vmc_step!(r, vmc_params, x -> ϕ(x, params))    
            E = EL(r, ϕ, params) 
            E_ave += E
        end
        
        add_sample!(Estat, E_ave/nn)
        for p in vmc_params
            adjust_step!(p)
        end

        if Estat.finished            
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)
            if verbose
                @printf("<E> =  %.5f +/- %.5f\n", Eb_ave, Eb_std)
            end
            # when to exit VMC
            if (Eb_std<tol && Nb>5)
                return Eb_ave, Eb_std
            end
        end        
    end
end


# Correlated local energy 
@inline function EL_correlated(walker ::Vector{Walker}, φ ::Function, params_new ::Params)
    E = 0.0        
    r = 0.0
    
    for w ∈ walker
        if w.Ψ<1e-12
            error("walker wave function w.Ψ is practically zero!")
        end
        wratio = φ(w.R, params_new)^2 / w.Ψ^2
        E += EL(w.R, φ, params_new) * wratio
        r += wratio
    end
    E = E/r
end

function generate_walkers(Nw ::Int64, Nwx ::Int64, ϕ ::Function, params:: Params, vmc_params :: Vector{VMC_Params})
    println("generating $Nw walkers, allocating $Nwx walkers")
    # electrons and protons
    r = [MMatrix{dim,N}(rand(dim,N) .- 0.5) for i in 1:Nwx]

    
    pdata = nothing
    walker = [Walker(;R=r[i], Ψ = ϕ(r[i], params), E=0.0, alive=0, age=0.0, weight=1.0, pure_data=pdata) for i in 1:Nwx]

    
    println("thermalizing walkers")
    E = 0.0
    alives = Set(1:Nw)  # set of alive walkers
    for iw in 1:Nw
        for i in 1:100
            vmc_step!(walker[iw].R, vmc_params, x -> ϕ(x, params))
            if i%10 == 0
                for p in vmc_params                    
                    adjust_step!(p)
                end                
            end
        end
        walker[iw].Ψ = ϕ(walker[iw].R, params)
        walker[iw].E = EL(walker[iw].R, ϕ, params)
        walker[iw].alive = 1
        E += walker[iw].E

    end
        
    for p in vmc_params
        @show p.step
    end
    E /= Nw
    println("<E> = $E")
    return walker, E, alives
end



# helpers for linesearch
# julia has no built-in pointers to elements of a structure, so create closures
# (getter and setter) to modify parameters

function getter(p::Params, field::Symbol)
    return () -> getfield(p, field)
end
function setter!(p::Params, field::Symbol)
    return (new_value) -> setfield!(p, field, new_value)
end


function linesearch(walker ::Vector{Walker}, ϕ ::Function, p ::Params, param_to_optimize ::Symbol)

    tol = 1e-8
    max_iter=10
    get_p = getter(p, param_to_optimize)
    set_p! = setter!(p, param_to_optimize)
    
    Nw = length(walker)
    Eold = EL_correlated(walker, ϕ, p)
    Δx = 0.1
    dir = 1
    for iter in 1:max_iter
        x = get_p()
        x += dir*Δx         
        if x < 0.1
            println("bad parameter")
            break
        end
        set_p!(x)
        Enew = EL_correlated(walker, ϕ, p)
        if Enew>Eold
            x -= dir*Δx
            set_p!(x)
            dir *= -1
            Δx *= 0.5
            if Δx<tol
                break
            end
        else
            Eold = Enew
            println("$Nw, $Enew, $p")
            open(outfile_opt,"a") do f
                println(f,"$Nw, $Enew, $p")
            end
        end            
    end
end

function lnG(xp ::MMatrix{dim,N,Float64},
             x ::MMatrix{dim,N,Float64},
             F ::MMatrix{dim,N,Float64},
             ϕ ::Function) where {dim, N}
    #
    # only parts not symmetric in x'<->x
    # G(x'<-x,τ) = exp(-(x'-x-DτF(x))^2/(4Dτ)
    #     
    get_F!(x, F, ϕ)
    lnG = 0.0
    Nhalf = Int(N/2)
    # electrons
    for k in 1:Nhalf
        dx = xp[:,k] - x[:,k] - D[1]*τ*F[:,k]
        lnG += sum(-dx.^2)/(4*D[1]*τ) 
    end
    # protons
    for k in Nhalf+1:N
        dx = xp[:,k] - x[:,k] - D[2]*τ*F[:,k]
        lnG += sum(-dx.^2)/(4*D[2]*τ) 
    end
    lnG
end

@noinline function get_F!(x ::MMatrix{dim,N,Float64}, F ::MMatrix{dim,N,Float64}, ϕ ::Function) where {dim, N}
    # for possible thread safety issues:
    #config = ForwardDiff.GradientConfig(ϕ, x) # a thread-local config makes gradient thread safe
    #F .= 2*ForwardDiff.gradient(ϕ, x, config)/ϕ(x)
    # 
    # Input F needs to be different for each thread!
    F .= 2*ForwardDiff.gradient(ϕ, x)/ϕ(x)
end

# walker update 1st order
@noinline function diffusion_drift_step!(x ::MMatrix{dim,N,Float64}, F ::MMatrix{dim,N,Float64}, ϕ ::Function) where {dim, N}  
    #
    rng = rngs[Threads.threadid()] 
    # diffusion(τ)+drift(τ)
    η = reshape(rand(rng, Normal(0,1), dim*N), (dim,N))
    get_F!(x, F, ϕ)
    Nhalf = Int(N/2)
    # electrons
    for k in 1:Nhalf
        @. x[:,k] += sqrt(2*D[1]*τ)*η[:,k] + D[1]*τ*F[:,k]
    end
    # protons
    for k in Nhalf+1:N
        @. x[:,k] += sqrt(2*D[2]*τ)*η[:,k] + D[2]*τ*F[:,k]
    end
    
        
end


function Ψ_i_AD(ϕ ::Function, R ::MMatrix{dim, N, Float64}, params ::Params) where {dim, N}  
    # φ_T(r ::MMatrix, p ::Params)
    para_concat = [params.ζ, params.b, params.c, params.d]   
    gradients = ForwardDiff.gradient(para_concat) do x
        params = Params(x[1], x[2], x[3], x[4])
        return ϕ(R, params)
    end
    return gradients
end

function EL_i_AD(ϕ ::Function, R ::MMatrix{dim, N, Float64}, params ::Params) where {dim, N}  
    # EL(r ::MMatrix, ϕ ::Function, params ::Params)
    para_concat = [params.ζ, params.b, params.c, params.d]
    gradients = ForwardDiff.gradient(para_concat) do x
        params = Params(x[1], x[2], x[3], x[4])
        return EL(R, ϕ, params)
    end
    return gradients
end

function optimization!(ϕ ::Function, vmc_params ::Vector{VMC_Params}; Nw_start=1000, Nw_stop=100000,
                       mode = "linear")
    # starting values
    ζ = 1.
    b = 0.5
    c = 1.
    d = 8.0
    params = Params(ζ, b, c, d)
    #
    # optimize parameters based on correlated local energy
    if mode=="linear"
        npara = 4
        Δα = @MVector zeros(npara)
        nonlin = [k for k in 1:npara]
        Nw = Nw_start
        while true # Nw loop 
            a_opt = 100.0
            

            Ψ_is  = zeros(npara, Nw)
            EL_is = zeros(npara, Nw)
            ELs   = zeros(Nw)

            for i in 1:5
                walker, Eold, alives = generate_walkers(Nw, Nw, ϕ, params, vmc_params)
                
                for k ∈ 1:10                    
                    for iw ∈ 1:Nw
                        R = walker[iw].R
                        Ψ_is[:,iw]  = Ψ_i_AD(ϕ, R, params) 
                        EL_is[:,iw] = EL_i_AD(ϕ, R, params) 
                        ELs[iw] =  EL(R, ϕ, params)
                        walker[iw].E = ELs[iw]
                        walker[iw].Ψ = ϕ(R, params)                    
                    end
                    Eold = sum(ELs)/Nw
                    
                    while true
                        solve_Δα!(walker, Ψ_is, EL_is, ELs, a_opt, npara, nonlin, Δα)
                        if maximum(abs.(Δα))<1e-1
                            break
                        end
                        a_opt *= 10.0
                    end
                        
                    params = update_params(params, Δα)
                    #
                    if params.ζ<0 || params.b<0 || params.c<0 || params.d<0
                        println("negative wf params")
                        params = update_params(params, -Δα)
                        a_opt *= 10.0
                    else
                        EL_corr_ave = EL_correlated(walker, ϕ, params)
                        if EL_corr_ave <= Eold
                            a_opt *= 0.5
                        else
                            params = update_params(params, -Δα)
                            a_opt *= 10.0                    
                        end
                        @printf("%5d a_opt %15.5e <EL> %15.7f <EL>_corr_ave %15.7f ",k,a_opt, Eold,  EL_corr_ave)
                        @printf("params %20.15f  %20.15f %20.15f %20.15f\n",params.ζ, params.b, params.c, params.d)
                        open(outfile_opt,"a") do f
                            println(f,"$Nw, $Eold, $params")
                        end
                    end
                    if a_opt>1e10 
                        break
                    end
                end
            end
            Nw = Int(floor(Nw*1.5))
            if Nw>Nw_stop
                break
            end
        end            
    else 
        params_to_optimize = [:ζ, :b, :c, :d]        
        Eold = 0.0
        Nw = Nw_start
        for i in 1:100
            walker, Eold, alives = generate_walkers(Nw, Nw, ϕ, params, vmc_params)
            for param_opt in params_to_optimize
                linesearch(walker, ϕ , params, param_opt)
            end
            Nw = floor(Int64,1.5*Nw)
            if Nw>Nw_stop
                break
            end
        end
    end
    return params
end

function main()

    vmc_params = [VMC_Params(0,0,1.0), VMC_Params(0,0,0.2)] # electrons, protons
    
    ϕ = φ_T  # trial wave function

    #
    # Linear or line optimization (choose mode="linear" or not)
    #
    if optimize==1
        params = optimization!(ϕ, vmc_params ; Nw_start=10000, Nw_stop=500000, mode="linear")        
        println("Equilibrium parameters $params")
        println("running VMC...")

        r = MMatrix{dim,N}(randn(dim,N))
        E, err = run_VMC!(r, ϕ, vmc_params, params ; tol=5e-4, verbose=true)
        @printf(" equilibrium <E> =  %.5f +/- %.5f\n", E, err)
        exit()
    end
    # Fix wave function parameters for VMC and DMC
    # linear optimization:
    #params = Params(1.2826976978333526, 0.3950542918239669, 1.3992489588778474, 11.057831033737072)
    # VMC : <E> =  -1.14627 +/- 0.00036
    #params = Params(1.2831514610054862, 0.39639318835767146, 1.399214103945839, 10.532560313210828)
    # VMC: <E> =  -1.14573 +/- 0.00029
    #
    # line search optimization:
    params = Params(1.28613, 0.39824, 1.3988, 10.5125)
    # VMC: <E> =  -1.14551 +/- 0.00012
    #params = Params(1.2796875, 0.396875, 1.403125, 10.80625)
    # VMC: <E> =  -1.14486 +/- 0.00035

    # params = Params(1.28242, 0.39805, 1.3885, 11.24375)
    #
    # =====
    #  VMC 
    # =====
    
    println("VMC with pre-optimized values")
    Nw = 10
    walker, E_ave, alives = generate_walkers(Nw, Nw, ϕ, params, vmc_params)    

    @show alives
    VMC_blocksize = 100 #10000 # fixed blocksize
    VMC_stat = init_stat_collection!(VMC_blocksize)
    # VMC end criterium:
    # make sure the most infrequent measurement gets a few blocks
    max_blocks = 50*maximum(values(measure_freq)) 
    @show max_blocks
    
    Estat = init_stat(1, VMC_blocksize) 
    E_VMC = 0.0
    ivmc = 1
    while true
        Req = 0.0
        E = 0.0
        @inbounds for iw in 1:Nw
            vmc_step!(walker[iw].R, vmc_params, x -> ϕ(x, params))
            E += EL(walker[iw].R, ϕ, params)  # measure energy
        end
        add_sample!(Estat, E/Nw)
        do_measurements(walker, ivmc, VMC_stat)

        # output finished measurements
        for measure in measurement_list
            if VMC_stat[measure].finished
                ave, std, _, Nb = get_stats(VMC_stat[measure])
                if data_size[measure]==1
                    @printf("VMC <%s> = %15.10f +/- %15.10f\n", string(measure), ave, std)
                elseif measure==:dens_profile
                    save_density_profile("VMC", ave, std)                    
                end
                # usually done in add_sample, but now samples are added infrequently
                VMC_stat[measure].finished = false
            end
        end

        if Estat.finished 
            E_VMC, E_std, _, Nb  = get_stats(Estat)
            @printf("VMC   <E> = %15.10f +/- %15.10f\n", E_VMC, E_std)
            println(" block $Nb of $max_blocks done")
            println("="^80)
            # stop VMC:
            if Nb==max_blocks
                break
            end
        end
        ivmc += 1
        
    end
    
    # =====
    #  DMC
    # =====
    println(("="^80)*"\n   DMC\n"*("=")^80)
    # drift; separate F for each thread
    Fs = [MMatrix{dim, N, Float64}(undef) for i in 1:Threads.nthreads()]

    # check drift
    vmc_params = [VMC_Params(0,0,1.0), VMC_Params(0,0,0.2)]
    walker, E_ave, alives = generate_walkers(10, 10, ϕ, params, vmc_params)

    @Threads.threads for iw = 1:10
        F = Fs[Threads.threadid()]
        get_F!(walker[iw].R, F , x->ϕ(x, params))
        # use x->F, a function that returns drift matrix F
        num_check_∇S(walker[iw].R, x->ϕ(x, params), x->F)
    end
    println("drift check passed")
    
    Nw = Nw_target
    walker, E_ave, alives = generate_walkers(Nw, Nwx, ϕ, params, vmc_params)
    ET = E_ave # trial energy, to be updated
    nE = 1 
    println("ET = $ET")
    
    
    idmc = -Ntherm
        
    copies = @MVector zeros(Int64,Nwx)
    
    Ntry = Threads.Atomic{Int64}(0)
    Nacc = Threads.Atomic{Int64}(0)
    
    # init measurements
    Estat = init_stat(1, blocksize)
    DMC_stat = init_stat_collection!(blocksize)

    while true
        Threads.@threads for iw in collect(alives) # vector from a set for threads
            #@inbounds for iw in alives
            rng = rngs[Threads.threadid()]
            F = Fs[Threads.threadid()]
            
            ELold = walker[iw].E
            Rnew = walker[iw].R # alias
            Rold  = copy(Rnew) # copy
            φTold = walker[iw].Ψ
            #
            diffusion_drift_step!(Rnew, F, x->ϕ(x, params))
            #
            Threads.atomic_add!(Ntry, 1)
            φTnew = ϕ(Rnew, params)

            Wnew = lnG(Rold, Rnew, F, x->ϕ(x, params)) + 2*log(φTnew)
            Wold = lnG(Rnew, Rold, F, x->ϕ(x, params)) + 2*log(φTold)
            # Metropolis
            accept = true
            if Wnew<Wold 
                if rand(rng)>exp(Wnew-Wold)
                    accept = false
                end
            end
            
            if accept
                Threads.atomic_add!(Nacc, 1)
                ELnew = EL(Rnew, ϕ, params)    
                walker[iw].E = ELnew
                walker[iw].Ψ = φTnew
            else
                ELnew = ELold
                Rnew .= Rold
            end
            weight  = exp(-τ*(0.5*(ELold+ELnew)-ET))
            copies[iw]  = floor(Int64, weight + rand(rng))
            walker[iw].weight = weight
        end
        branch!(walker, copies, alives)
        Nw = length(alives)
        E   = 0.0
        @inbounds for iw in alives
            E += walker[iw].E                                     
        end
        E_ave += E/Nw
        nE += 1
        idmc +=1
        
        # update trial energy
        ET = E_ave/nE + κ*log(Nw_target/Nw)

        # measurements
        if idmc>=0
            add_sample!(Estat, E/Nw)
            do_measurements(walker, idmc, DMC_stat)
        end
        open(Eallfile,"a") do f
            println(f, E/Nw," ",E_ave/nE," " ,ET," ",Nw)
        end
        # collect atomic values
        total_Ntry = Ntry.value
        total_Nacc = Nacc.value
        if idmc<0
            @printf("DMC %10d  τ = %.5e <E> = %.10f E = %.10f ET = %.10f  %6d Walkers  acceptance = %.5f \n",
                    idmc, τ, E_ave/nE,  E/Nw, ET, Nw, 100.0*total_Nacc/total_Ntry)
        end
        # block data
        # ==========
        # screen and file output

        for measure in measurement_list
            if DMC_stat[measure].finished
                ave, std, _, Nb = get_stats(DMC_stat[measure])
                ave_VMC, std_VMC, _, _ = get_stats(VMC_stat[measure])      
                if data_size[measure]==1                                  
                    @printf("    VMC <%s> = %15.10f +/- %.10f\n", string(measure), ave_VMC, std_VMC)
                    @printf("    DMC <%s> = %15.10f +/- %.10f\n", string(measure), ave, std)
                    @printf("Extrap. <%s> = %15.10f +/- %.10f\n", string(measure), 2*ave-ave_VMC, sqrt((2*std)^2+(std_VMC)^2))
                elseif measure==:dens_profile
                    save_density_profile("DMC", ave, std)
                    extrap = @.  2*ave-ave_VMC
                    sigma_extrap = @. sqrt((2*std)^2+(std_VMC)^2)
                    save_density_profile("Extrap", extrap, sigma_extrap)
                end
                DMC_stat[measure].finished = false
            end
        end
        
        
        if Estat.finished            
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)

            @printf("DMC %10d τ = %.5e  <E> = %.10f +/- %.10f E = %.10f ET = %.10f \
                %6d Walkers \n", idmc, τ, Eb_ave, Eb_std,  E/Nw,  ET, Nw)
            
            open(Efile,"a") do f
                println(f,τ," ",Eb_ave," ",Eb_std," ",Nw)
            end
            
            if Nb>10 && Eb_std < accuracy_goal
                println("reached accuracy goal $accuracy_goal")
                println("result: τ  <E>  error")
                @printf("%g ",τ)
                open(Efile,"a") do f
                    println(f,τ," ",Eb_ave," ",Eb_std)
                end
                output_MCresult(Eb_ave, Eb_std)
                break
            end
            GC.gc() # manual carbage collection
        end
        
        # 
        if idmc==0            
            println("THERMALIZATION ENDS")
            #Estat = init_stat(1, blocksize) # only needed if measuring during thermalization
            E_ave = E_ave/nE
            nE = 1
        end
    end
end

main()


