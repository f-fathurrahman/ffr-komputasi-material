
#
# H2 molecule
# Trial wave function from Kolos and Roothaan, Rev. Mod. Phys. 32 205 (1960) 
#
# For example, run as
# (i) optimization of parameters for proton-proton distance Req=1.2:
#  julia H2_Kolos_Roothaan.jl optimize=1 Req=1.2
# (ii) DMC, using time step 0.001 and proton-proton distance Req=1.2 
#  julia H2_Kolos_Roothaan.jl tau=0.001 Req=1.2 
#
#  default Req=1.4 



using Printf
using StaticArrays
import Distributions: Normal
using Statistics
using ForwardDiff
using ForwardDiff: gradient, derivative, Dual, value
#using ReverseDiff: hessian
using LinearAlgebra: diag, tr
using Random


@inline norm(x) = sqrt(sum(abs2,x))

push!(LOAD_PATH,".")
using Common: VMC_Params, Walker, branch!
using LinearOptimization
using QMC_Statistics
using Utilities: argparse, num_check_EL, num_check_∇S, output_MCresult
using VMCstep: adjust_step!, vmc_step_H2! as vmc_step!

const dim=3
const N = 4
const D = 1/2


# Threads and RNG
println("Using $(Threads.nthreads()) threads")
# seed rng's using time in milliseconds
base_number = floor(Int, time() * 1000) # floor needed to avoid InexactError
#println("FIXED SEED TESTING")
#base_number = floor(Int, 12345 * 1000) # floor needed to avoid InexactError
#base_number = 12345
rngs = [MersenneTwister(base_number + i) for i in 1:Threads.nthreads()]


# Wave function parameters
mutable struct Params{T<:Real}
    α ::T    
    a ::Vector{T}
    b ::Vector{T}
    c ::Vector{T}
end

# Print parameters in copy-paste format
function print_params(p::Params)
    print("params = Params(")
    @printf("%8.6f,\n", p.α)
    
    l = length(p.a)
    for s in [p.a, p.b, p.c]
        print("[")
        for i in 1:l
            @printf("%8.6f", s[i])
            if i < l
                print(", ")  
            end
        end        
        print("]")
        if s in [p.a, p.b]
            println(",")
        end
    end
    println(")")
end


function update_params!(p ::Params, Δp ::MVector)
    p.α += Δp[1]
    l = length(p.a)
    offs = 2
    p.a .+= Δp[offs:offs+l-1]
    p.b .+= Δp[offs+l:offs+2l-1]
    p.c .+= Δp[offs+2l:end]
end

# Command line parameters
# either 
#   optimize=1
# or
#   tau= (give imaginary time step) 
# 
possible_args=["optimize","tau", "Req"]
arg_dict = argparse(possible_args)
optimize = get(arg_dict, "optimize", 0)
Req_in =  get(arg_dict, "Req", 1.4)
const Req = Req_in
@show Req

if optimize==1
    
    const outfile_opt = "H2_Kolos_Roothaan_opt.dat"
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
    const blocksize = floor(Int64,1/(5τ)) :: Int64 # just an approx
    const Ntherm = floor(Int64,2/(5τ))

   
    @show blocksize
    @show Ntherm
    
    const Efile = "E.H2_Kolos_Roothaan_"*string(τ)*"_"*@sprintf("%.2f", Req)
    @show Efile
    rm(Efile, force=true)
    const Eallfile = "E.tmp_"*string(τ)*"_"*@sprintf("%.2f", Req)
    rm(Eallfile, force=true)
  
end

# measurements
# ============

const measurement_list = [] # [:Vpot, :dens_profile]
# sizes of measured data
const data_size = Dict(:Vpot=>1, :dens_profile=>301)
# frequencies of measurements, expensive measurements infrequent
const measure_freq = Dict(:Vpot=>5, :dens_profile=>10)

dens_grid = range(-3.0, 4.0, length=data_size[:dens_profile])

# Actual measurements


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
#
# Kolos and Roothaan Rev. Mod. Phys. 32 205 (1960)
# 
# Eq.(5)
function u_Kolos_Rothaan(ξ ::T, η ::U, pi ::Int, qi ::Int, α ::V) where {T,U,V}
    # q_i must be even
    u = ξ^pi * η^qi * exp(-α*ξ) 
end

function v_Kolos_Rothaan(r12 ::T, μ ::Int) where T
    v = r12^μ
end

# Table I
const p_i = [0, 1, 2, 0, 1]
const q_i = [0, 0, 0, 2, 2]
const μ_i = [0,1,2,3,4] 

function φ_Kolos_Rothaan(r ::MMatrix, p ::Params) 
    α, a, b, c = p.α, p.a, p.b, p.c 
    r1, r2, r3, r4 = r[:,1], r[:,2], r[:,3], r[:,4]
    r13 = norm(r1-r3)
    r14 = norm(r1-r4)
    r23 = norm(r2-r3)
    r24 = norm(r2-r4)
    r12 = norm(r1-r2)
    r34 = norm(r3-r4)
    # elliptic coordinates
    R = Req

    ξ1 = (r13+r14)/R
    ξ2 = (r23+r24)/R
    η1 = (r13-r14)/R
    η2 = (r23-r24)/R

    u = u_Kolos_Rothaan
    v = v_Kolos_Rothaan
    
    # sums in Eq.(4)
    ϕ1 = 0.0
    ϕ2 = 0.0
    ψ1 = 0.0
    ψ2 = 0.0
    for i in 1:length(a)
        ϕ1 += a[i] * u(ξ1, η1, p_i[i], q_i[i], α)
        ϕ2 += a[i] * u(ξ2, η2, p_i[i], q_i[i], α)        
        ψ1 += b[i] * u(ξ1, η1, p_i[i], q_i[i], α)
        ψ2 += b[i] * u(ξ2, η2, p_i[i], q_i[i], α)
    end
    χ12 = 0.0
    for i in 1:length(c)
        χ12 += c[i] * v(r12, μ_i[i])  
    end
    wf = (ϕ1*ψ2 + ψ1*ϕ2) * χ12 # Eq.(3)

end



function φ_T(r ::MMatrix, p ::Params)    
    wf = φ_Kolos_Rothaan(r, p)
end

# local kinetic energy TL =-D  ∇^2 ϕ(r) / ϕ(r)
function TL(r, ϕ ::Function)   
    sum_∇2ϕ = 0.0
    dim, N = size(r)
    Nhalf = Int(N/2)
    for i in 1:Nhalf # electron TL 
        for k in 1:dim
            # FIXME: this sets also proton r's to Dual, which works but is unnecessary
            function f(t)
                r_copy = similar(r, typeof(t)) 
                r_copy .= r  
                r_copy[k, i] = t                
                return ϕ(r_copy)      
            end
            sum_∇2ϕ += derivative(t -> derivative(f, t), r[k, i])
        end
    end
    return -D*sum_∇2ϕ/ϕ(r) 
    
end

# potential energy
@inline function V(r ::MMatrix) 
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
function EL(r, ϕ ::Function, params ::Params)
    return TL(r , x -> ϕ(x, params)) + V(r)    
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

function generate_walkers(Nw ::Int64, Nwx ::Int64, ϕ ::Function, p:: Params, vmc_params :: VMC_Params)
    println("generating $Nw walkers, allocating $Nwx walkers")
    # electrons and protons
    r = [MMatrix{dim,N}(rand(dim,N) .- 0.5) for i in 1:Nwx]
    for iw in 1:Nwx    
        r[iw][:, 3] = 0.0
        r[iw][:, 4] = 0.0
        r[iw][1, 4] = Req
    end
    
    pdata = nothing
    walker = [Walker(;R=r[i], Ψ = ϕ(r[i], p), E=0.0, alive=0, age=0.0, weight=1.0, pure_data=pdata) for i in 1:Nwx]

    
    println("thermalizing walkers")
    E = 0.0
    alives = Set(1:Nw)  # set of alive walkers
    for iw in 1:Nw
        R = walker[iw].R
        for i in 1:300
            vmc_step!(R, vmc_params, x -> ϕ(x, p))
            if i%10 == 0
                adjust_step!(vmc_params)
            end
        end
        walker[iw].Ψ = ϕ(R, p)
        walker[iw].E = EL(R, ϕ, p)
        walker[iw].alive = 1
        E += walker[iw].E
    end
        
    @show vmc_params.step
    E /= Nw
    println("<E> = $E")
    return walker, E, alives
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
        dx = xp[:,k] - x[:,k] - D*τ*F[:,k]
        lnG += sum(-dx.^2)/(4*D*τ) 
    end
    lnG
end

@noinline function get_F!(x ::MMatrix{dim,N,Float64}, F ::MMatrix{dim,N,Float64}, ϕ ::Function) where {dim, N}
    # for possible thread safety issues:
    #config = ForwardDiff.GradientConfig(ϕ, x) # a thread-local config makes gradient thread safe
    #F .= 2*ForwardDiff.gradient(ϕ, x, config)/ϕ(x)
    # 
    # Input F needs to be different for each thread!
    F .= 2*gradient(ϕ, x)/ϕ(x)
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
        @. x[:,k] += sqrt(2*D*τ)*η[:,k] + D*τ*F[:,k]
    end
end


function Ψ_i_AD(ϕ ::Function, R ::MMatrix{dim, N, Float64}, p ::Params) where {dim, N}  
    # φ_T(r ::MMatrix, p ::Params)
    para_concat = [p.α]
    l = length(p.a)
    for i in 1:l
        append!(para_concat, p.a[i])
    end
    for i in 1:l
        append!(para_concat, p.b[i])
    end
    for i in 1:l
        append!(para_concat, p.c[i])
    end
    gradients = gradient(para_concat) do x
        p = Params(x[1], x[2:2+l-1], x[2+l:2+2l-1], x[2+2l:2+3l-1]) 
        return ϕ(R, p)
    end
    return gradients
end

function EL_i_AD(ϕ ::Function, R ::MMatrix{dim, N, Float64}, p ::Params) where {dim, N}  
    # EL(r ::MMatrix, ϕ ::Function, p ::Params)
    para_concat = [p.α]
    l = length(p.a)
    for i in 1:l
        append!(para_concat, p.a[i])
    end
    for i in 1:l
        append!(para_concat, p.b[i])
    end
    for i in 1:l
        append!(para_concat, p.c[i])
    end
    gradients = gradient(para_concat) do x
        p = Params(x[1], x[2:2+l-1], x[2+l:2+2l-1], x[2+2l:2+3l-1]) 
        return EL(R, ϕ, p)
    end
    return gradients
end

function optimization!(ϕ ::Function, p ::Params, vmc_params ::VMC_Params; Nw_start=1000, Nw_stop=100000,
                       mode = "linear")
    
    # optimize parameters based on correlated local energy
    if mode=="linear"
        npara = 1+3*length(p.a)
        Δα = @MVector zeros(npara)
        nonlin = [k for k in 1:npara]
        Nw = Nw_start
        while true # Nw loop 
            a_opt = 100.0
            
            Ψ_is  = zeros(npara, Nw)
            EL_is = zeros(npara, Nw)
            ELs   = zeros(Nw)

            for i in 1:1
                walker, Eold, alives = generate_walkers(Nw, Nw, ϕ, p, vmc_params)
                Eorig = Eold
                for k ∈ 1:10                  
                    for iw ∈ 1:Nw
                        R = walker[iw].R
                        Ψ_is[:,iw]  = Ψ_i_AD(ϕ, R, p) 
                        EL_is[:,iw] = EL_i_AD(ϕ, R, p) 
                        ELs[iw] =  EL(R, ϕ, p)
                        walker[iw].E = ELs[iw]
                        walker[iw].Ψ = ϕ(R, p)                    
                    end
                    Eold = sum(ELs)/Nw
                    if (abs(Eorig-Eold)>0.1)
                        break
                    end
                    while true
                        solve_Δα!(walker, Ψ_is, EL_is, ELs, a_opt, npara, nonlin, Δα)
                        @show Δα
                        Δα[2:end] .*= 0.01 # slow down
                        #if maximum(abs.(Δα))<1e-2                            
                        break
                        #end
                        #a_opt *= 2.0
                    end                        
                    update_params!(p, Δα)
                    #
                    #if p.ζ<0 || p.b<0 || p.c<0 || p.d<0 || p.η[1]<0|| p.η[2]<0
                    #    println("negative wf p")
                    #    update_params!(p, -Δα)
                    #    a_opt *= 10.0
                    if false
                        continue
                    else
                        EL_corr_ave = EL_correlated(walker, ϕ, p)
                        if EL_corr_ave <= Eold + 1/Nw^2
                            #Eold = EL_corr_ave # try this
                            a_opt *= 0.1
                        else
                            update_params!(p, -Δα)
                            a_opt *= 10.0                    
                        end
                        @printf("%5d a_opt %15.5e <EL> %15.7f <EL>_corr_ave %15.7f \n",k,a_opt, Eold,  EL_corr_ave)
                        print_params(p)    
                        
                        open(outfile_opt,"a") do f
                            println(f,"$Nw, $Eold, $p")
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
        for i in 1:length(p.η)
            push!(params_to_optimize, Symbol("η_$i"))
        end
        @show params_to_optimize 

        Eold = 0.0
        Nw = Nw_start
        for i in 1:100
            walker, Eold, alives = generate_walkers(Nw, Nw, ϕ, p, vmc_params)
            for param_opt in params_to_optimize
                linesearch(walker, ϕ , p, param_opt, 0.1/Nw)
            end
            Nw = floor(Int64,1.5*Nw)
            if Nw>Nw_stop
                break
            end
        end
    end
    return p
end

# pre-optimized parameters
function preoptimized_params(Req ::Float64)
    params_dict = Dict()

    params_dict[0.02] = Params(0.033116,
                               [1.019979, -0.033840, 0.000969, 0.012674, 0.000061],
                               [1.024255, 0.024836, 0.000165, 0.017531, -0.000144],
                               [0.921213, 0.440283, -0.095365, 0.014785, -0.000932])
    
    params_dict[0.04] = Params(0.052677,
                               [1.019973, -0.034027, 0.001017, 0.012680, 0.000044],
                               [1.024256, 0.024770, 0.000319, 0.017537, -0.000130],
                               [0.921211, 0.440288, -0.095369, 0.014758, -0.001007])

    params_dict[0.06] = Params(0.073433,
                               [1.019957, -0.034511, 0.000956, 0.012699, 0.000024],
                               [1.024257, 0.024688, 0.000931, 0.017547, -0.000105],
                               [0.921211, 0.440289, -0.095369, 0.014736, -0.001117])

    params_dict[0.1] = Params(0.099498,
                              [1.019919, -0.035615, 0.000884, 0.012724, -0.000126],
                              [1.024260, 0.024495, 0.001815, 0.017560, -0.000097],
                              [0.921211, 0.440290, -0.095353, 0.014781, -0.001030])    

    params_dict[0.2] = Params(0.184471,
                              [1.019842, -0.037671, 0.000762, 0.012987, -0.000243],
                              [1.024268, 0.023857, 0.003508, 0.017712, -0.000056],
                              [0.921166, 0.440411, -0.095202, 0.014963, -0.000967])
    
    params_dict[0.3] = Params(0.265276,
                              [1.019807, -0.038535, 0.000569, 0.013227, -0.000278],
                              [1.024259, 0.023713, 0.004233, 0.018246, 0.000090],
                              [0.921090, 0.440771, -0.094301, 0.014725, -0.000830])
    
    params_dict[0.75] = Params(0.583518,
                               [1.005905, 0.045598, 0.014644, 0.087538, 0.004132],
                               [1.004134, -0.085978, 0.002466, 0.075262, -0.005520],
                               [0.920041, 0.444642, -0.063139, 0.006811, -0.000266])
    
    params_dict[1.0] = Params(0.741498,
                              [0.995800, 0.035919, 0.023542, 0.152719, 0.002589],
                              [0.994977, -0.094219, 0.002024, 0.115571, -0.006148],
                              [0.918901, 0.448724, -0.047122, 0.003267, -0.000014])

    params_dict[1.4] = Params(0.984346,
                              [0.988673, -0.031269, 0.047094, 0.233969, 0.015059],
                              [0.985987, -0.072300, -0.003928, 0.225481, -0.016072],
                              [0.932543, 0.456358, -0.017962, -0.002725, 0.000404])

    params_dict[2.0] = Params(1.316357,
                              [0.950864, -0.209694, 0.082444, 0.321109, 0.026997],
                              [0.918931, 0.039663, -0.029757, 0.462699, -0.030198],
                              [0.927152, 0.438568, 0.029403, -0.008576, 0.000589])
    
    params_dict[2.5] = Params(1.560596,
                              [0.891703, -0.260430, 0.079224, 0.448178, 0.015413],
                              [0.839395, 0.019713, -0.023777, 0.608471, -0.042931],
                              [0.927800, 0.434986, 0.056675, -0.007749, 0.000347])    

    params_dict[3.0] = Params(1.898421,
                              [0.860375, -0.415552, 0.118183, 0.469255, 0.061443],
                              [0.641672, 0.301747, -0.090116, 0.836668, -0.064068],
                              [0.943542, 0.391003, 0.111150, -0.010599, 0.000547])

    params_dict[3.5] = Params(2.247789,
                              [0.753122, -0.392769, 0.103167, 0.682585, 0.067262],
                              [0.762227, 0.195645, -0.056098, 0.829498, -0.025148],
                              [0.964784, 0.332928, 0.168619, -0.011085, 0.000844])

    params_dict[4.0] = Params(2.614946,
                              [0.728304, -0.415089, 0.099145, 0.696144, 0.073873],
                              [0.753529, 0.196254, -0.053496, 0.838196, 0.004889],
                              [0.968992, 0.313345, 0.181942, -0.006673, 0.001167])

    params_dict[4.5] = Params(3.042140,
                              [0.698120, -0.434638, 0.096318, 0.714104, 0.084693],
                              [0.734570, 0.186720, -0.058793, 0.856594, 0.012824],
                              [0.970091, 0.301311, 0.196054, -0.001891, 0.001868])

    params_dict[5.0] =  Params(3.390702,
                               [0.684363, -0.457477, 0.092887, 0.712712, 0.092207],
                               [0.713892, 0.168192, -0.054564, 0.878063, 0.008580],
                               [0.978958, 0.247877, 0.225757, -0.002920, 0.003052])

    params_dict[6.0] = Params(4.225819,
                              [0.647531, -0.504296, 0.097400, 0.715735, 0.086342],
                              [0.680396, 0.138826, -0.053832, 0.908295, 0.045177],
                              [0.977221, 0.194095, 0.279980, -0.009131, 0.006733])

    params_dict[7.0] = Params(4.932698,
                              [0.609589, -0.562428, 0.118582, 0.703877, 0.075340],
                              [0.552416, 0.136187, -0.058993, 0.992374, 0.044828],
                              [0.936373, -0.217069, 0.397916, -0.042617, 0.009941])
    
    params_dict[10.0] =Params(7.102220,
                              [0.541324, -0.619918, 0.148102, 0.710260, -0.014007],
                              [0.500120, -0.062579, -0.020479, 1.027099, 0.089721],
                              [0.648708, -0.688157, 0.447878, -0.088061, 0.007896])


    
    # Find closes parameter set to current Req
    allkeys = collect(keys(params_dict))
    key = allkeys[argmin(abs.(allkeys .- Req))]
    println("Using pre-optimized parameters at Req=$key")  
    return params_dict[key]
    
end

function main()

    vmc_params = VMC_Params(0,0,1.0)
    
    ϕ = φ_T  # trial wave function
    
    params = preoptimized_params(Req)
    print_params(params)
    
    #
    # Linear optimization
    #
    if optimize==1       
        #        
        params = optimization!(ϕ, params, vmc_params ; Nw_start=1000, Nw_stop=500000)        
        println("Equilibrium parameters")
        print_params(params)        
        exit()
    end

    
    #
    # =====
    #  VMC 
    # =====
    
    println("VMC with pre-optimized values")
    Nw = 10
    walker, E_ave, alives = generate_walkers(Nw, Nw, ϕ, params, vmc_params)    

    @show alives
    VMC_blocksize = 1000 # fixed blocksize
    VMC_stat = init_stat_collection!(VMC_blocksize)
    # VMC end criterium:
    # make sure the most infrequent measurement gets a few blocks
    max_blocks = 10*maximum(values(measure_freq)) 
    @show max_blocks
    
    Estat = init_stat(1, VMC_blocksize) 
    E_VMC = 0.0
    ivmc = 1

    # save plot of electron cloud
    if false
        Nsave = 1000
        Rfile = "R_"*@sprintf("%.2f", Req)
        println("Saving electron coordinates to file $Rfile")
        open(Rfile,"w") do f
            for i in 1:Nsave
                @inbounds for iw in 1:Nw
                    vmc_step!(walker[iw].R, vmc_params, x -> ϕ(x, params))
                end     
                @inbounds for iw in 1:Nw
                    @printf(f,"%f8.5  %f8.5   %f8.5\n", walker[iw].R[:,1]...)
                    @printf(f,"%f8.5  %f8.5   %f8.5\n", walker[iw].R[:,2]...)                
                end
            end
            Nsave += 1
        end
        exit()
    end
    
    while true
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
    vmc_params = VMC_Params(0,0,1.0)
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


