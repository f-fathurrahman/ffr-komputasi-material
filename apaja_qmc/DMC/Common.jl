# Common to many QMC programs

module Common
using StaticArrays
using ForwardDiff: Dual
using ForwardDiff
using ReverseDiff
using DelimitedFiles: readdlm
using CurveFit

import Base: one, zero

# slightly faster norm 
@inline norm(x) = sqrt(sum(abs2,x))

export Walker, init_walkers, VMC_Params, spinfactor
export branch!

export AtomicOrbital, LinearCombination, JastrowParameters
export dim, D
export V

export parse_atom_data, parse_basis_set_data, parse_Jastrow_data, parse_STO, parse_coeffients_data
export WfParams
export update_wf_params,  norm_params!
export set_linear_combination_of_orbitals

const dim = 3    # space dimension
const D = 0.5    # hbar^2/2m in a.u. 


struct AtomData
    atom::String
    N::Int
    Nup::Int
    Ndo::Int
    Z::Int
    Eexact::Float64
end

# use two Number types, because AD derives exponents and coefficients separately
# ForwardDiff uses Dual type, ReverseDiff uses TrackedArray type
struct AtomicOrbital{A<:AbstractString, T<:Number, Z<:Number}
    atom         ::A
    orbital_type ::A
    n            ::Int
    l            ::Int
    exponents    ::AbstractVector{T}
    coefficients ::AbstractVector{Z}
    spin         ::A
end

struct LinearCombination{T<:Number}
    orbitals     ::Vector{Any}
    coeffs       ::AbstractVector{T}
end

struct Molecule{T<:Number,Z<:Number}
    atoms        ::Vector{String}          # Atomic species, e.g., ["H", "H", "O"]
    positions    ::Vector{SVector{3, T}}   # Atom positions; 3D 
    molecular_orbitals ::Vector{LinearCombination{Z}}  # MO's (LCAO)
end

function Ψ_mol(m::Molecule, R::MMatrix{dim,N,Float64}) where {dim, N}
    Ψ = 0.0
    for mo in m.molecular_orbitals
        for (ao, coeff) in zip(mo.orbitals, mo.coeffs)
            # Add contributions from atomic orbitals
            Ψ += coeff * DD(ao, R)
        end
    end
    return Ψ
end

function DD(ao::AtomicOrbital, R::MMatrix{dim,N,Float64}) where {dim, N}
    DD = 0.0    
    for (i,aos) in enumerate(lc.orbitals)
        ups = aos[lc_ups[i]]
        dos = aos[lc_dos[i]]
        DD += lc.coeffs[i] * D_up(R, ups) * D_do(R, dos)
    end
end

struct JastrowParameters{T<:Number}
    param ::AbstractVector{T}
end

struct WfParams
    lc ::LinearCombination
    jpara ::JastrowParameters
end

# define copy of wave function parameters
Base.copy(p::WfParams) = WfParams(p.lc, p.jpara)

# normalize c's so that the 1s is always 1.0
# (also lc.coeffs may start with 1)
function norm_params!(wf_params::WfParams, basis_set ::String, Z ::Int64)
    
    sc = 0.0
    for aos in wf_params.lc.orbitals        
        for ao in aos
            if ao.orbital_type == "1s"
                sc = ao.coefficients[1]
                break
            end
        end
    end
    
    orbitals = []    
    for aos in wf_params.lc.orbitals
        aolist = []
        for ao in aos                
            new_coefficients = ao.coefficients ./ sc
            new_ao = AtomicOrbital(ao.atom, ao.orbital_type, ao.n, ao.l, ao.exponents, new_coefficients, ao.spin)
            push!(aolist, new_ao)
        end
        push!(orbitals, aolist)        
    end
    
    # normalize lc.coeffs
    scc =  wf_params.lc.coeffs[1]
    wf_params.lc.coeffs ./= scc
    
    # keep 2p coefficients same    
    n_2p = 0
    s_2p = 0.0
    for (aos, coef) in zip(wf_params.lc.orbitals, wf_params.lc.coeffs)
        for ao in aos
            if occursin("2p", ao.orbital_type)                
                s_2p += coef
                n_2p += 1
            end
        end
    end
    # average
    s_2p /= n_2p
    # modify lc.coeffs that contain 2p orbitals
    i = 0
    for aos in wf_params.lc.orbitals
        i += 1
        for ao in aos            
            if occursin("2p", ao.orbital_type)
                wf_params.lc.coeffs[i] = s_2p
            end
        end
    end
   
    
    lc = LinearCombination(orbitals, wf_params.lc.coeffs)
    jpara = wf_params.jpara
    return  WfParams(lc, jpara)

end

function one(::Type{AtomicOrbital})
    return AtomicOrbital("", "", 0,0,ones(Float64, 0), ones(Float64, 0), 0)
end


function zero(::Type{AtomicOrbital})
    return AtomicOrbital("", "", zeros(Float64, 0), zeros(Float64, 0), 0)
end

import Base: +

function +(ao::AtomicOrbital, b::AtomicOrbital) 
        AtomicOrbital(ao.atom, ao.orbital_type, ao.n, ao.l, ao.exponents .+ b.exponents, ao.coefficients .+ b.coefficients, ao.spin)
end

# In case of ReverseDiff, define + between TrackedArray and StaticArray types
function +(x::ReverseDiff.TrackedArray, y::StaticArray)
    return x .+ y  # Element-wise broadcasting
end

function +(x::StaticArray, y::ReverseDiff.TrackedArray)
    return x .+ y  # Element-wise broadcasting
end

# define WfParams + Δα
# Notice that Δα is in a specific order
#  
function update_wf_params(a::WfParams, x::MVector{T}) where T
    lc = a.lc    
    pos = 1
    len = length(lc.coeffs)
    new_coeffs = lc.coeffs .+ x[pos:pos+len-1]
    pos += len
    
    nl = []
    new_params = Dict()
    for aos in lc.orbitals
        for ao in aos            
            if [ao.n, ao.l] ∉ nl
                push!(nl,[ao.n, ao.l])
                len = length(ao.exponents)
                # don't let exponents change more than, say, 10 %
                
                new_exponents = ao.exponents .+ x[pos:pos+len-1]
                # don't let exponents change more than, say, 10 %
                for iexp ∈ 1:len
                    dd = new_exponents[iexp]-ao.exponents[iexp]
                    if abs(dd)/abs(ao.exponents[iexp])>0.10
                        new_exponents[iexp] = ao.exponents[iexp] + 0.10*dd
                    end
                end
                
                # check that exponents are positive
                if !all(>=(0),new_exponents)
                    println("Common update_wf_params! : negative orbital exponent")  
                    @show ao.exponents, new_exponents
                    return a
                end
                pos += len
                new_coefficients = ao.coefficients .+ x[pos:pos+len-1]
                pos += len
                new_params[ao.n,ao.l] = [new_exponents,new_coefficients]
            end
        end
    end
    
    orbitals = []
    for aos in lc.orbitals
        new_aolist = AtomicOrbital[]                
        for ao in aos            
            new_ao = AtomicOrbital(ao.atom, ao.orbital_type, ao.n, ao.l, new_params[ao.n,ao.l][1], new_params[ao.n,ao.l][2], ao.spin)
            push!(new_aolist, new_ao)            
        end
        push!(orbitals, new_aolist)
    end
    lc = LinearCombination(orbitals, new_coeffs)

    param = copy(a.jpara.param)  # don't set param=a.jpara.param, you'll chance input!
   
    for i ∈ 1:4
        param[i] +=  x[end-4+i]
    end
    # check that Jastrow β1 and β2 are positive
    if param[2]<0.0 ||  param[3]<0.0
        println("Common update_wf_params! : negative Jastrow β1 or β2")  
        #jpara = a.jpara.param
        return a
    end
        
    jpara = JastrowParameters(param)

    new_a = WfParams(lc, jpara)
    return new_a
end


mutable struct VMC_Params
    Ntry :: Int64
    Naccept :: Int64
    step :: Float64
end 


mutable struct Walker
    R     :: MMatrix
    alive :: Int64
    Ψ     :: Float64    
    E     :: Float64
    age   :: Float64
    weight :: Float64
    pure_data ::Union{Nothing, T} where T
    measurements :: Dict{Symbol, Any} # other measured quantities
    # constructor
    Walker(; R, alive, Ψ, E, age, weight, pure_data = nothing) =
        new(R, alive, Ψ,  E, age, weight, pure_data, Dict{Symbol, Any}())
end


# 
# Initialize walkers
# ==================
# Simplified version
# 
function init_walkers(Nw ::Int64, dim ::Int64, N ::Int64) ::Array{Walker,1}

    println("Init $Nw walkers  - note: most values set to zero")
    MMzeros = @MMatrix zeros(dim,N)
    # take copy of MMzeros, else the R's reference the same memory location
    walker = [Walker(R = copy(MMzeros), alive=0, Ψ=0.0, E=0.0, age=0.0,
                     weight=1.0) for i in 1:Nw]

    walker
end

# Initialization with trial wave function
function init_walkers(Nw ::Int64,
                      Nwx ::Int64,
                      dim ::Int64,
                      N ::Int64;
                      φ_T  ::Function = dummy,
                      scale ::Float64 = 1.0,                      
                      ) ::Array{Walker,1}
    #
    println("generating $Nw walkers, storage for $Nwx walkers")
    MMzeros = @MMatrix zeros(dim,N)
    # take copy of MMzeros, else the R's reference the same memory location
    walker = [Walker(R = copy(MMzeros), alive=0, Ψ=0.0, E=0.0, age=0.0, weight=1.0) for i in 1:Nwx]

    # set Nw walkers
    for w ∈ walker[1:Nw]
        w.alive = 1
        while true
            R = @MMatrix rand(dim,N)
            R .*= scale # scale a bit
            Ψ = φ_T(R)
            if Ψ > 0
                w.R = copy(R)  
                w.Ψ = copy(Ψ)
                break
            end
        end        
    end
    walker
end

    

# Walker to Walker copy
@inline function copy!(dest::Walker, src::Walker)
    if dest.R === src.R
        error("dest and src R's reference the same memory location.")
    end
    dest.R .= src.R    # broadcast
    dest.alive = src.alive
    dest.Ψ = src.Ψ    
    dest.E = src.E
    dest.age = src.age
    dest.weight = src.weight
    if src.pure_data !== nothing
        dest.pure_data = deepcopy(src.pure_data) 
    else
        dest.pure_data = nothing  
    end
    dest.measurements =  deepcopy(src.measurements)
end



@inline function V(R ::MMatrix, N ::Int64, Z ::Int64) ::Float64
    V = 0.0 ::Float64
    for k in 1:N
        rk  = norm(R[:,k])
        V += -Z/rk
        for i in k+1:N
            rki = norm(R[:,k]-R[:,i]) 
            V += 1/rki
        end
    end
    V
end




#
# Set up CSF's (Configuration State Functions)
# - done manually in the atom_data setup files, nothing is automated
#
function set_linear_combination_of_orbitals(atom_data ::AtomData,
                                            aos ::Vector{AtomicOrbital},
                                            cnls ::Vector
                                            )
       

    println("="^80)
    println("="^80)
    println("Calculation will use these atomic orbitals:") 
    
    atom = atom_data.atom
    
    # collect requested atomic orbitals listed in cnls
    aolist = []
    orbital_dict = Dict( (1,0,0)=>"1s", (2,0,0)=>"2s", (2,1,-1)=>"2px", (2,1,0)=>"2py" ,(2,1,1)=>"2pz")                     

    coeffs = Float64[]
    m = -1
    for i in 1:2:length(cnls)
        coeff = cnls[i]
        append!(coeffs, coeff)
        nls = cnls[i+1]
        spin="up"
        
        aolist_ = AtomicOrbital[]

        for j in 1:2:length(nls)
            n = nls[j]
            l = nls[j+1]
            # find correct ao from known aos
            found = false
            for (i,ao) in enumerate(aos)
                if ao.n==n && ao.l==l                    
                    orbital_type = orbital_dict[(n,l,0)]
                    if n==2 && l==1
                        orbital_type = orbital_dict[(n,l,m)]
                    end
                    
                    ao_ = AtomicOrbital(ao.atom,
                                        orbital_type,
                                        ao.n,
                                        ao.l,
                                        ao.exponents,
                                        ao.coefficients,
                                        spin)
                    push!(aolist_,ao_)                                            
                    found = true
                end
            end            
            if !found
                error("Atomic Orbital n l = $n $l is not among known aos.")
            end
            # next spin
            spin = (spin == "up") ? "down" : "up"
            if n==2 && l==1 && spin=="up"
                # next m
                if m==-1
                    m = 0
                elseif m==0
                    m = 1
                else
                    m = -1
                end
            end

        end
        push!(aolist,aolist_)
        
    end
    
    # only different n,l pairs are orbitals with different parameters
    nl = []
    npara = length(coeffs) # coefficients are first parameters
    nparanonlin = 0
   
    for (coeff,aos) in zip(coeffs,aolist)
        println("="^80)
        @show coeff
        for ao in aos
            @show ao
            if [ao.n, ao.l] ∉ nl
                push!(nl, [ao.n, ao.l])
                npara += 2*length(ao.exponents)
                nparanonlin += length(ao.exponents)
            end
        end
    end
    println("="^80)

    lc = LinearCombination(
        aolist, # list of atomic orbitals 
        coeffs # Coefficients        
    )

    nonlin = Int64[]
    pos = length(lc.coeffs) # first params are lc.coeffs
    nl = []
    for aos in lc.orbitals
        for ao in aos        
            if [ao.n, ao.l] ∉ nl
                push!(nl,[ao.n, ao.l])
                len = length(ao.exponents)
                append!(nonlin, [i for i in pos+1:pos+len])
                pos += 2*len # skip coefficients
            end
        end
    end
    # list Jastrow α, β1 and β2 as nonlinear parameters ; not yet counted in npara
    append!(nonlin, [npara+1, npara+2, npara+3])
            
    @show nonlin
    lc_ups = []
    lc_dos = []
    for aos in lc.orbitals
        ups = Int[]
        dos = Int[]
        for (i,ao) in enumerate(aos)
            if ao.spin =="up"
                push!(ups,i)
            else
                push!(dos,i)
            end
        end
        push!(lc_ups,ups)
        push!(lc_dos,dos)
    end
    return lc, npara, nonlin, lc_ups, lc_dos
   
end


# ====================

# spin factor
function spinfactor(Nup, Ndo)
    N = Nup+Ndo
    s = zeros(Int64,N,N)
    for i ∈ 1:N
        for j ∈ 1:N
            s[i,j] = spins(i, j, Nup, Ndo)
        end
    end
    s
end

function spins(i ::Int64, j ::Int64, Nup ::Int64, Ndo ::Int64) ::Int64    
    if (i<= Nup && j>Nup) || (i> Nup && j<=Nup)
        # anti-parallel spins                
        fac = 2 ::Int64    
    else
        # parallel spins
        fac = 4 ::Int64    
    end
    fac
end

function print_basis_set(basis_set_data)
    for ao in basis_set_data
        println("Orbital Type: ", ao.orbital_type)
        try
            println("Exponents: ", [ζ.value for ζ in ao.exponents])
            println("Coefficients: ", [c.value for c in ao.coefficients])
        catch
            println("Exponents: ", ao.exponents)
            println("Coefficients: ", ao.coefficients)
        end
    end
end


function dummy(x ::Any)
    return 1.0
end






# Common DMC functions
# ====================


function branch!(walker ::Vector{Walker}, copies ::MVector, alives ::Set{Int64})
    Nwx = length(walker)
    # list initially dead walkers that can we can branch to
    deads = findall(iw -> iw ∉ alives, 1:Nwx)
    idead = 1
    orig_alives = copy(alives) # modify alives in loop
    @inbounds for iw in orig_alives
        if copies[iw]==1
            continue        # one copy is already there
        elseif copies[iw]==0
            # walker dies
            walker[iw].alive = 0
            delete!(alives, iw)
            continue
        end
        # copies[iw]>1
        # copy the walker to empty slots 
        @inbounds for inew in 1:copies[iw]-1   # with copies=3 this is 1,2
            # use free slot at idead
            if idead > length(deads)
                @show length(alives)
                error("No free walkers available; too large τ ?")
            end
            iw2 = deads[idead]
            copy!(walker[iw2], walker[iw]) # walker to walker copy
            push!(alives, iw2)             
            idead += 1 # take next free slot 
        end
    end
    return nothing
end


# ================
# Parsers
# =================





function parse_basis_set_data(atom::String, basis_set ::String)
    STOfile = "atom_data/$(atom)_$(basis_set)"
    if isfile(STOfile)
        basis_set_data = parse_STO(STOfile)
    else
        error("File not found: $STOfile")
    end
    #
    #    # Fall back to STO-3g data and fit them to STO
    #    println("Warning: Failed to read data from $STOfile, falling back to STO-3g")
    #    error("not tested")
    #    STO3gfile = "atom_data/"*atom*"_STO3g"
    #    basis_set_data = parse_turbomole_sto3g(STO3gfile)
    #    basis_set_data, Spara, nonlin = fit_STO3g_to_STO(basis_set_data)
    #    Norb = 0 # not needed
    #end
    basis_set_data
end

function parse_Jastrow_data(atom::String, basis_set ::String) ::JastrowParameters
    filename = "atom_data/"*atom*"_"*basis_set
    println("Reading Jastrow wf parameters from file $filename")
    α, β1, β2, α12 = 0.0, 0.0, 0.0, 0.0
    open(filename, "r") do file
        while !eof(file)
            line = strip(readline(file))
            # Remove comments
            line = split(line, "#")[1]
            line = strip(line)
            if startswith(line, "α =")
                α = parse(Float64, strip(split(line, "=")[2]))
            elseif startswith(line, "β1 =")
                β1 = parse(Float64, strip(split(line, "=")[2]))
            elseif startswith(line, "β2 =")
                β2 = parse(Float64, strip(split(line, "=")[2]))
            elseif startswith(line, "α12 =")
                α12 = parse(Float64, strip(split(line, "=")[2]))
            end
        end
    end
    para = [α, β1, β2, α12]
    return JastrowParameters(para)
end


function parse_coeffients_data(atom::String, basis_set ::String)
    filename = "atom_data/"*atom*"_"*basis_set
    println("Reading linear combination coefficients from file $filename")
    cnls=[]
    open(filename, "r") do file
        while !eof(file)
            line = strip(readline(file))
            # Remove comments
            line = split(line, "#")[1]
            line = strip(line)
            if startswith(line, "coeff =")                
                cc = parse(Float64, strip(split(line, "=")[2]))
                append!(cnls,cc)
            end
            nl = Int[]
            if startswith(line, "orbitals =")
                parts = split(line)               
                for p in parts[3:end]
                    snl = string(parse(Int,p))
                    n = parse(Int, snl[1])
                    l = parse(Int, snl[2])
                    append!(nl,n,l)
                end
                append!(cnls,[nl])  
            end
        end
    end
    # output looks like
    # cnls = Any[1.0, [1, 0, 1, 0, 2, 0, 2, 0], -0.4464125065, [1, 0, 1, 0, 2, 1, 2, 1],
    #            -0.4464125065, [1, 0, 1, 0, 2, 1, 2, 1], -0.4464125065, [1, 0, 1, 0, 2, 1, 2, 1]]

    return cnls
end



# Function to parse the data from the file
function parse_atom_data(atom::String)
    N = 0
    Nup = 0
    Ndo = 0
    Z = 0
    Eexact = 0.0
    

    filename = "atom_data/"*atom
    println("reading atom data from file $filename")
    open(filename, "r") do file
        while !eof(file)
            line = strip(readline(file))
            # Remove comments
            line = split(line, "#")[1]
            line = strip(line)
            
            if isempty(line)
                continue
            end
            if startswith(line, "N =")
                N = parse(Int, strip(split(line, "=")[2]))
            elseif startswith(line, "Nup =")
                Nup = parse(Int, strip(split(line, "=")[2]))
            elseif startswith(line, "Ndo =")
                Ndo = parse(Int, strip(split(line, "=")[2]))
            elseif startswith(line, "Z =")
                Z = parse(Int, strip(split(line, "=")[2]))
            elseif startswith(line, "Eexact =")
                Eexact = parse(Float64, strip(split(line, "=")[2]))               
            end
        end
    end
    
    return AtomData(atom, N, Nup, Ndo, Z, Eexact)
end

# Function to parse the atomic orbital data from STO file  
function parse_STO(filename ::String)
    l_dict  = Dict("s"=>0, "p"=>1, "d"=>2, "e"=>3, "f"=>4, "g"=>5, "h"=>6, "i"=>7, "j"=>8)
    fill_order = [1,2,2,3,3,4,3,4,5,4,5,6,4,5,6]
    
    
    println("Reading atomic orbital data from file $filename")
    atom = ""  
    orbitals = AtomicOrbital[]
    i = 0
    println("="^80)
    println("Known orbitals and parameters:") 
    open(filename, "r") do file
        while !eof(file)
            line = strip(readline(file))
            
            if startswith(line, "#")
                # Skip comment lines
                continue
            elseif startswith(line, "*")
                # End of a section
                continue
            elseif isempty(atom)
                # First non-comment line is the atom
                atom = string(split(line)[1]) #just to avoid being a substring
            else
                # Parsing the orbital data
                parts = split(line)
                
                if length(parts) == 2 && parse(Int, parts[1]) > 0
                    # Start of a new orbital section
                    num_primitives = parse(Int, parts[1])
                    orbital_type = parts[2]
                    
                    
                    exponents = Float64[]
                    coefficients = Float64[]
                    
                    for _ in 1:num_primitives
                        line = strip(readline(file))
                        values = split(line)
                        push!(exponents, parse(Float64, values[1]))
                        push!(coefficients, parse(Float64, values[2]))
                    end
                    i += 1
                    n = fill_order[i]
                    l = l_dict[orbital_type]
                    s = string("$n"*orbital_type)
                    orbital =  AtomicOrbital(atom, s, n, l, exponents, coefficients, "spin_not_set")
                    @show orbital
                    push!(orbitals, orbital)
                end
            end 
        end
    end
    return orbitals 
end


# Function to parse the atomic orbital data from the Turbomole file
function parse_turbomole_sto3g(STO3gfile::String)
    orbitals = []
    filename = STO3gfile
    println("Reading atomic orbital data from file $filename")
    atom = ""
    
    open(filename, "r") do file
        while !eof(file)
            line = strip(readline(file))
            
            if startswith(line, "#")
                # Skip comment lines
                continue
            elseif startswith(line, "*")
                # End of a section
                continue
            elseif isempty(atom)
                # First non-comment line is the atom
                atom = split(line)[1]
            else
                # Parsing the orbital data
                parts = split(line)
                
                if length(parts) == 2 && parse(Int, parts[1]) > 0
                    # Start of a new orbital section
                    num_primitives = parse(Int, parts[1])
                    orbital_type = parts[2]
                    
                    exponents = Float64[]
                    coefficients = Float64[]
                    
                    for _ in 1:num_primitives
                        line = strip(readline(file))
                        values = split(line)
                        push!(exponents, parse(Float64, replace(values[1], "D" => "E")))
                        push!(coefficients, parse(Float64, replace(values[2], "D" => "E")))
                    end
                    
                    orbital =  AtomicOrbital(atom, orbital_type, exponents, coefficients)
                    push!(orbitals, orbital)
                end
            end
        end
    end
    
    return orbitals
end
#=
# Example usage

orbitals = parse_turbomole_sto3g(filename)

# Print the parsed data
for orbital in orbitals
println("Atom: ", orbital.atom)
println("Orbital Type: ", orbital.orbital_type)
println("Exponents: ", orbital.exponents)
println("Coefficients: ", orbital.coefficients)
println()
end
=#


end

