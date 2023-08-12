# Binary bases commonly used in canonical fermionic and hardcore bosonic quantum lattice systems
"""
    BinaryBasis{I<:Unsigned}

Binary basis represented by an unsigned integer.
"""
struct BinaryBasis{I<:Unsigned}
    rep::I
    BinaryBasis(i::Integer) = (rep = Unsigned(i); new{typeof(rep)}(rep))
end
@inline Base.:(==)(basis₁::BinaryBasis, basis₂::BinaryBasis) = basis₁.rep == basis₂.rep
@inline Base.isequal(basis₁::BinaryBasis, basis₂::BinaryBasis) = isequal(basis₁.rep, basis₂.rep)
@inline Base.:<(basis₁::BinaryBasis, basis₂::BinaryBasis) = basis₁.rep < basis₂.rep
@inline Base.isless(basis₁::BinaryBasis, basis₂::BinaryBasis) = isless(basis₁.rep, basis₂.rep)
@inline Base.one(basis::BinaryBasis) = one(typeof(basis))
@inline Base.one(::Type{BinaryBasis{I}}) where {I<:Unsigned} = BinaryBasis(one(I))
@inline Base.zero(basis::BinaryBasis) = zero(typeof(basis))
@inline Base.zero(::Type{BinaryBasis{I}}) where {I<:Unsigned} = BinaryBasis(zero(I))
@inline Base.show(io::IO, basis::BinaryBasis) = @printf io "%s" string(basis.rep, base=2)
@inline Base.eltype(basis::BinaryBasis) = eltype(typeof(basis))
@inline Base.eltype(::Type{<:BinaryBasis}) = Int
@inline Base.IteratorSize(::Type{<:BinaryBasis}) = Base.SizeUnknown()






"""
    iterate(basis::BinaryBasis, state=nothing)

Iterate over the numbers of the occupied single-particle orbitals.
"""
function Base.iterate(basis::BinaryBasis, state=nothing)
    (pos, rep) = isnothing(state) ? (0, basis.rep) : (state[1], state[2])
    while rep>0
        pos += 1
        isodd(rep) && return (pos, (pos, rep÷2))
        rep ÷= 2
    end
    return nothing
end

"""
    one(basis::BinaryBasis, state::Integer) -> BinaryBasis

Get a new basis with the specified single-particle state occupied. 
"""
@inline Base.one(basis::BinaryBasis, state::Integer) = BinaryBasis(basis.rep | one(basis.rep)<<(state-1))

"""
    isone(basis::BinaryBasis, state::Integer) -> Bool

Judge whether the specified single-particle state is occupied for a basis.
"""
@inline Base.isone(basis::BinaryBasis, state::Integer) = (basis.rep & one(basis.rep)<<(state-1))>0

"""
    zero(basis::BinaryBasis, state::Integer) -> BinaryBasis

Get a new basis with the specified single-particle state unoccupied.
"""
@inline Base.zero(basis::BinaryBasis, state::Integer) = BinaryBasis(basis.rep & ~(one(basis.rep)<<(state-1)))

"""
    iszero(basis::BinaryBasis, state::Integer) -> Bool

Judge whether the specified single-particle state is unoccupied for a basis.
"""
@inline Base.iszero(basis::BinaryBasis, state::Integer) = !isone(basis, state)

"""
    count(basis::BinaryBasis) -> Int
    count(basis::BinaryBasis, start::Integer, stop::Integer) -> Int

Count the number of occupied single-particle states.
"""
@inline Base.count(basis::BinaryBasis) = count(basis, 1, ndigits(basis.rep, base=2))
@inline function Base.count(basis::BinaryBasis, start::Integer, stop::Integer)
    result = 0
    for i = start:stop
        isone(basis, i) && (result += 1)
    end
    return result
end

"""
    ⊗(basis₁::BinaryBasis, basis₂::BinaryBasis) -> BinaryBasis

Get the direct product of two binary bases.
"""
@inline ⊗(basis₁::BinaryBasis, basis₂::BinaryBasis) = BinaryBasis(basis₁.rep|basis₂.rep)

"""
    BinaryBasis(states; filter=index->true)
    BinaryBasis{I}(states; filter=index->true) where {I<:Unsigned}

Construct a binary basis with the given occupied orbitals.
"""
@inline BinaryBasis(states; filter=index->true) = BinaryBasis{typeof(Unsigned(first(states)))}(states; filter=filter)
function BinaryBasis{I}(states; filter=index->true) where {I<:Unsigned}
    rep, eye = zero(I), one(I)
    for (index, state) in enumerate(states)
        filter(index) && (rep += eye<<(state-1))
    end
    return BinaryBasis(rep)
end









"""
    abstract type Sector <: OperatorUnit

A sector of the Hilbert space which form the bases of an irreducible representation
of the Hamiltonian of a quantum lattice system.
"""
abstract type Sector <: OperatorUnit end







"""
    BinaryBasisRange{I<:Unsigned} <: VectorSpace{BinaryBasis{I}}

A continuous range of binary basis.
"""
struct BinaryBasisRange{I<:Unsigned} <: VectorSpace{BinaryBasis{I}}
    slice::UnitRange{I}
end
@inline Base.issorted(bbr::BinaryBasisRange) = true
@inline Base.length(bbr::BinaryBasisRange) = length(bbr.slice)
@inline Base.getindex(bbr::BinaryBasisRange, i::Integer) = BinaryBasis(bbr.slice[i])








"""
    BinaryBases{B<:BinaryBasis, T<:AbstractVector{B}} <: Sector

A set of binary bases.
"""
struct BinaryBases{B<:BinaryBasis, T<:AbstractVector{B}} <: Sector
    id::Vector{Tuple{B, Float64}}
    table::T
end
@inline Base.issorted(bs::BinaryBases) = true
@inline Base.length(bs::BinaryBases) = length(bs.table)
@inline Base.:(==)(bs₁::BinaryBases, bs₂::BinaryBases) = isequal(bs₁.id, bs₂.id)
@inline Base.isequal(bs₁::BinaryBases, bs₂::BinaryBases) = isequal(bs₁.id, bs₂.id)
@inline Base.getindex(bs::BinaryBases, i::Integer) = bs.table[i]
@inline Base.eltype(bs::BinaryBases) = eltype(typeof(bs))
@inline Base.eltype(::Type{<:BinaryBases{B}}) where {B<:BinaryBasis} = B
@inline Base.iterate(bs::BinaryBases, state=1) = state>length(bs) ? nothing : (bs.table[state], state+1)
function Base.repr(bs::BinaryBases)
    result = String[]
    for (rep, nparticle) in bs.id
        if isnan(nparticle)
            push!(result, @sprintf "2^%s" count(rep))
        else
            push!(result, @sprintf "C(%s, %s)" count(rep) Int(nparticle))
        end
    end
    return join(result, " ⊗ ")
end
function Base.show(io::IO, bs::BinaryBases)
    @printf io "%s:\n" repr(bs)
    for i = 1:length(bs)
        @printf io "  %s\n" bs[i]
    end
end
@inline Base.searchsortedfirst(b::BinaryBasis, bs::BinaryBases) = searchsortedfirst(bs.table, b)
@inline Base.searchsortedfirst(b::BinaryBasis, bs::BinaryBases{<:BinaryBasis, <:BinaryBasisRange}) = Int(b.rep+1)





"""
    ⊗(bs₁::BinaryBases, bs₂::BinaryBases) -> BinaryBases

Get the direct product of two sets of binary bases.
"""
function ⊗(bs₁::BinaryBases, bs₂::BinaryBases)
    @assert productable(bs₁, bs₂) "⊗ error: the input two sets of bases cannot be direct producted."
    table = Vector{promote_type(eltype(bs₁), eltype(bs₂))}(undef, length(bs₁)*length(bs₂))
    count = 1
    for (b₁, b₂) in product(bs₁, bs₂)
        table[count] = b₁⊗b₂
        count += 1
    end
    return BinaryBases(sort!([bs₁.id; bs₂.id]; by=first), sort!(table))
end

"""
    productable(bs₁::BinaryBases, bs₂::BinaryBases) -> Bool

Judge whether two sets of binary bases could be direct producted.
"""
function productable(bs₁::BinaryBases, bs₂::BinaryBases)
    for (irr₁, irr₂) in product(bs₁.id, bs₂.id)
        isequal(irr₁[1].rep & irr₂[1].rep, 0) || return false
    end
    return true
end

"""
    sumable(bs₁::BinaryBases, bs₂::BinaryBases) -> Bool

Judge whether two sets of binary bases could be direct summed.

Strictly speaking, two sets of binary bases could be direct summed if and only if they have no intersection. The time complexity to check the intersection is O(n log n), which costs a lot when the dimension of the binary bases is huge. It is also possible to judge whether they could be direct summed by close investigations on their ids, i.e. the single-particle states and occupation number. It turns out that this is a multi-variable pure integer linear programming problem. In the future, this function would be implemented based on this observation. At present, the direct summability should be handled by the users in priori.
"""
@inline sumable(bs₁::BinaryBases, bs₂::BinaryBases) = true

"""
    BinaryBases(nstate::Integer)
    BinaryBases(states)

Construct a set of binary bases that does not preserve the particle number conservation.
"""
@inline BinaryBases(nstate::Integer) = BinaryBases([(BinaryBasis(1:nstate), NaN)], BinaryBasisRange(UInt(0):UInt(2^nstate-1)))
@inline BinaryBases(states) = BinaryBases(Tuple(states))
function BinaryBases(states::NTuple{N, Integer}) where N
    states = NTuple{N, eltype(states)}(sort!(collect(states); rev=true))
    com = DuplicatePermutations{N}((false, true))
    table = Vector{BinaryBasis{typeof(Unsigned(first(states)))}}(undef, length(com))
    for (i, poses) in enumerate(com)
        table[i] = BinaryBasis(states; filter=index->poses[index])
    end
    return BinaryBases([(BinaryBasis(states), NaN)], table)
end

"""
    BinaryBases(nstate::Integer, nparticle::Integer)
    BinaryBases(states, nparticle::Integer)

Construct a set of binary bases that preserves the particle number conservation.
"""
@inline BinaryBases(nstate::Integer, nparticle::Integer) = BinaryBases(1:nstate, Val(nparticle))
@inline BinaryBases(states, nparticle::Integer) = BinaryBases(states, Val(nparticle))
function BinaryBases(states, ::Val{N}) where N
    com = Combinations{N}(sort!(collect(states); rev=true))
    I = typeof(Unsigned(first(states)))
    table = Vector{BinaryBasis{I}}(undef, length(com))
    for (i, poses) in enumerate(com)
        table[end+1-i] = BinaryBasis{I}(poses)
    end
    return BinaryBases([(BinaryBasis{I}(states), Float64(N))], table)
end






"""
    EDKind{K}

The kind of the exact diagonalization method applied to a quantum lattice system.
"""
struct EDKind{K} end
@inline EDKind(K::Symbol) = EDKind{K}()
@inline EDKind(::Type{T}) where {T<:Term} = error("EDKind error: not defined for $(kind(T)).")
@inline EDKind(::Type{<:FockTerm}) = EDKind(:FED)
#@inline EDKind(::Type{<:SpinTerm}) = EDKind(:SED) # ffr
@inline @generated function EDKind(::Type{TS}) where {TS<:Tuple{Vararg{Term}}}
    exprs = []
    for i = 1:fieldcount(TS)
        push!(exprs, :(typeof(EDKind(fieldtype(TS, $i)))))
    end
    return Expr(:call, Expr(:call, :reduce, :promote_type, Expr(:tuple, exprs...)))
end









"""
    ED{K<:EDKind, L<:AbstractLattice, G<:OperatorGenerator, M<:Image} <: Frontend

Exact diagonalization method of a quantum lattice system.
"""
struct ED{K<:EDKind, L<:AbstractLattice, G<:OperatorGenerator, M<:Image} <: Frontend
    lattice::L
    H::G
    Hₘ::M
    function ED{K}(lattice::AbstractLattice, H::OperatorGenerator, mr::EDMatrixRepresentation) where K
        Hₘ = mr(H)
        new{K, typeof(lattice), typeof(H), typeof(Hₘ)}(lattice, H, Hₘ)
    end
end
@inline kind(ed::ED) = kind(typeof(ed))
@inline kind(::Type{<:ED{K}}) where K = K()
@inline Base.valtype(::Type{<:ED{<:EDKind, <:AbstractLattice, G}}) where {G<:OperatorGenerator} = valtype(eltype(G))
@inline statistics(ed::ED) = statistics(typeof(ed))
@inline statistics(::Type{<:ED{<:EDKind, <:AbstractLattice, G}}) where {G<:OperatorGenerator} = statistics(eltype(eltype(G)))
@inline function update!(ed::ED; kwargs...)
    if length(kwargs)>0
        update!(ed.H; kwargs...)
        update!(ed.Hₘ, ed.H)
    end
    return ed
end
@inline Parameters(ed::ED) = Parameters(ed.H)



"""
    ED(lattice::AbstractLattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, targetspace::TargetSpace; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)

Construct the exact diagonalization method for a quantum lattice system.
"""
function ED(lattice::AbstractLattice, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, targetspace::TargetSpace; neighbors::Union{Nothing, Int, Neighbors}=nothing, boundary::Boundary=plain)
    k = EDKind(typeof(terms))
    isnothing(neighbors) && (neighbors = maximum(term->term.bondkind, terms))
    H = OperatorGenerator(terms, bonds(lattice, neighbors), hilbert; half=false, boundary=boundary)
    mr = EDMatrixRepresentation(targetspace, Table(hilbert, Metric(k, hilbert)))
    return ED{typeof(k)}(lattice, H, mr)
end



