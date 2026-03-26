module Utilities
using Printf
using StaticArrays

export output_MCresult, num_check_EL, num_check_∇S, metro
export show_vec
export argparse
export print_params
export dist1, dist2

@inline function dist2(R::MMatrix{dim, N, Float64}, i::Int64, j::Int64) where {dim, N}
    s::Float64 = 0.0
    @inbounds for k in 1:dim
        dx = R[k, i] - R[k, j]
        s += dx*dx
    end
    return sqrt(s)
end

@inline function dist1(R::MMatrix{dim, N, Float64}, i::Int64) where {dim, N}
    s::Float64 = 0.0
    @inbounds for k in 1:dim
        s += R[k, i]^2
    end
    return sqrt(s)
end

function print_params(para,txt=""::String)
    @printf("%10s ",txt)
    try
        for cc in para.lc.coeffs
            @printf("%8.5f ", cc)
        end
        nl = []
        for aos in para.lc.orbitals
            for ao in aos
                if [ao.n, ao.l] ∉ nl
                    push!(nl,[ao.n, ao.l]) 
                    for ζ in ao.exponents
                        @printf("%8.5f ",ζ)
                    end
                    for c in ao.coefficients
                        @printf("%8.5f ",c)
                    end
                end
            end
        end
        for p in para.jpara.param
            @printf("%8.5f ",p)
        end
        println()
    catch        
        try            
            for p in para
                @printf("%8.5f ",p)
            end
            println()
        catch
            println("can't use print_params")
        end
    end
end



# =====================================
# ArgParse would do this, but I hate to load yet another package

function argparse(known_args ::Vector{String})
    if length(ARGS)==0; return ; end

    got_args = Dict() 
    for arg in ARGS
        if occursin("=",arg)
            parts = split(arg, "=", keepempty=false)
            var = parts[1]
            if var ∉ known_args
                @show(known_args)
                println("Error reading command line arguments, variable $var is unknown")
                exit(1)
            end
            try
                value = parts[2]
                #@show(value)
                if occursin(".",value)
                    value = parse(Float64,value)
                else
                    try
                        value = parse(Int64,value)
                    catch
                        # assume it's a string
                        value = String(value)
                    end
                end
                #println("command line: ",var,"=",value)
                got_args[var]=value
            catch
                println("Error reading command line arguments, use var=val")
                exit(1)
            end
        end
    end
    return got_args
end
# ======================================



# result output as
# value +/- error
#  and
# value(error)
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
function output_MCresult(stream, value, error)
    if isapprox(error,0.0)
        acc =  floor(Int64, log10(1/1e-10)+2)        
    else
        acc = floor(Int64, log10(1/error)+2)
    end
    fmt = Printf.Format("%."*"$(acc)f"*" +/- "*"%."*"$(acc)f \n")
    Printf.format(stream, fmt, value, error)
end



# Metropolis question for ln-weights
@noinline function metro(Wold ::Float64, Wnew ::Float64)
    if Wnew>Wold
        return true
    end
    
    if rand()<exp(Wnew-Wold)
        return true
    else
        return false
    end
end

function show_vec(vec ::T) where T<:AbstractVector{Float64}    
    map(x->@printf("%10.6f ",x), vec)
    println()
end


# EL check
function num_check_EL(R ::MMatrix, EL_in, Ψ ::Function, V ::Function)
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

# ∇S check
function num_check_∇S(R ::MMatrix, Ψ ::Function, drift ::Function)
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






# ADAM algorithm
# https://arxiv.org/pdf/1412.6980.pdf
function adam!(θ ::Vector, acc ::Float64)
    n = length(θ)
    α = 0.1  # Stepsize
    β1 = 0.9 # Exponential decay rates for momentum estimates
    β2 = 0.999    
    ϵ = 1e-8    
    m = zeros(n) 
    v = zeros(n) 
    t = 0 # timestep
    while true
        t += 1
        g = grad
        m  =  β1*m + (1-β1)*g
        v  = β2*v + (1-β2)*(g.*g)
        mhat = m/(1-β1^t)
        vhat = v/(1-β2^t)
        upd = -α*mhat/(√vhat + ϵ)
        @. θ += update
        if maximum(abs.(update))<acc
            break
        end
    end
end


end
