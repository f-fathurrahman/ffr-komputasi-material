using Base.Iterators: product


## Constraint
const wildcard = Symbol("*")



"""
    DataType(T::DataType) -> DataType
    DataType(T::UnionAll) -> DataType

Get the DataType.
"""
@inline Base.DataType(T::DataType) = T
@inline Base.DataType(T::UnionAll) = DataType(T.body)




"""
    parameterorder(::Type{T}, name::Symbol) where T -> Int

For a type `T`, get the order of one of its type parameters.
"""
@inline parameterorder(::Type{T}, name::Symbol) where T = _order(Val(name), parameternames(T)|>Val)
@inline @generated _order(::Val{name}, ::Val{names}) where {name, names} = findfirst(isequal(name), names)::Int




"""
    parametertype(::Type{T}, name::Symbol) where T
    parametertype(::Type{T}, i::Integer) where T

For a type `T`, get the type of one of its type parameters.
"""
@inline parametertype(::Type{T}, name::Symbol) where T = parametertype(T, parameterorder(T, name))
@inline parametertype(::Type{T}, i::Integer) where T = _parametertype(T, Val(i))
@inline @generated function _parametertype(::Type{T}, ::Val{i}) where {T, i}
    result = DataType(T).parameters[i]
    return isa(result, TypeVar) ? result.ub : result
end





"""
    getcontent(m, i::Integer)
    getcontent(m, name::Symbol)
    getcontent(m, ::Val{name}) where name

Get the value of the predefined content of `m`. 
"""
@inline getcontent(m, i::Integer) = getcontent(m, contentname(typeof(m), i))
@inline getcontent(m, name::Symbol) = getcontent(m, Val(name))
@inline getcontent(m, ::Val{name}) where name = getfield(m, name)




"""
    tile(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations) -> Matrix{<:Number}

Tile a supercluster by translations of the input cluster.

Basically, the final supercluster is composed of several parts, each of which is a translation of the original cluster, with the translation vectors specified by `vectors` and each set of the translation indices contained in `translations`. When translation vectors are empty, a copy of the original cluster will be returned.
"""
function tile(cluster::AbstractMatrix{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}, translations)
    length(vectors)==0 && return copy(cluster)
    length(translations)>0 && @assert length(vectors)==length(first(translations)) "tile error: mismatched shape of input vectors and translations."
    datatype = promote_type(eltype(cluster), eltype(eltype(vectors)), Float)
    supercluster = zeros(datatype, size(cluster, 1), size(cluster, 2)*length(translations))
    disp = zeros(datatype, size(cluster, 1))
    for (i, translation) in enumerate(translations)
        for i = 1:length(disp)
            disp[i] = zero(datatype)
            for j = 1:length(vectors)
                disp[i] += vectors[j][i] * translation[j]
            end
        end
        for j = 1:size(cluster, 2)
            col = (i-1)*size(cluster, 2) + j
            for row = 1:size(cluster, 1)
                supercluster[row, col] = cluster[row, j] + disp[row]
            end
        end
    end
    return supercluster
end




"""
    AbstractLattice{N, D<:Number, M}

Abstract type of a unitcell-described lattice.

It should have the following contents:
- `name::Symbol`: the name of the lattice
- `coordinates::Matrix{D}`: the coordinates of the lattice
- `vectors::SVector{M, SVector{N, D}}`: the translation vectors of the lattice
"""
abstract type AbstractLattice{N, D<:Number, M} end
@inline contentnames(::Type{<:AbstractLattice}) = (:name, :coordinates, :vectors)
@inline Base.:(==)(lattice₁::AbstractLattice, lattice₂::AbstractLattice) = ==(efficientoperations, lattice₁, lattice₂)
@inline Base.isequal(lattice₁::AbstractLattice, lattice₂::AbstractLattice) = isequal(efficientoperations, lattice₁, lattice₂)
@inline Base.eltype(lattice::AbstractLattice) = eltype(typeof(lattice))
@inline Base.eltype(::Type{<:AbstractLattice{N, D}}) where {N, D<:Number} = SVector{N, D}
@inline Base.iterate(lattice::AbstractLattice, state=1) = state>length(lattice) ? nothing : (lattice[state], state+1)


function Base.show(io::IO, lattice::AbstractLattice)
    @printf io "%s(%s)\n" lattice|>typeof|>nameof getcontent(lattice, :name)
    len = length(lattice)
    if len > 0
        @printf io "  with %s %s:\n" len len==1 ? "point" : "points"
        for i = 1:len
            @printf io "    %s\n" lattice[i]
        end
    end
    len = length(getcontent(lattice, :vectors))
    if len > 0
        @printf io "  with %s translation %s:\n" len len==1 ? "vector" : "vectors"
        for i = 1:len
            @printf io "    %s\n" getcontent(lattice, :vectors)[i]
        end
    end
end



"""
    dimension(lattice::AbstractLattice) -> Int
    dimension(::Type{<:AbstractLattice{N}}) where N -> Int

Get the space dimension of the lattice.
"""
@inline dimension(lattice::AbstractLattice) = dimension(typeof(lattice))
@inline dimension(::Type{<:AbstractLattice{N}}) where N = N

"""
    dtype(lattice::AbstractLattice)
    dtype(::Type{<:AbstractLattice{N, D} where N}) where {D<:Number}

Get the data type of the coordinates of a lattice.
"""
@inline dtype(lattice::AbstractLattice) = dtype(typeof(lattice))
@inline dtype(::Type{<:AbstractLattice{N, D} where N}) where {D<:Number} = D

"""
    length(lattice::AbstractLattice) -> Int

Get the number of points contained in a lattice.
"""
@inline Base.length(lattice::AbstractLattice) = size(getcontent(lattice, :coordinates))[2]

"""
    getindex(lattice::AbstractLattice, i::Integer) -> SVector

Get the ith coordinate.
"""
@inline Base.getindex(lattice::AbstractLattice, i::Integer) = SVector{dimension(lattice), dtype(lattice)}(ntuple(j->lattice.coordinates[j, i], Val(dimension(lattice))))

"""
    reciprocals(lattice::AbstractLattice) -> Vector{<:SVector}

Get the reciprocal translation vectors of the dual lattice.
"""
@inline reciprocals(lattice::AbstractLattice) = reciprocals(getcontent(lattice, :vectors))







"""
    Lattice{N, D<:Number, M} <: AbstractLattice{N, D, M}

Simplest lattice.

A simplest lattice can be constructed from its coordinates and translation vectors.
"""
struct Lattice{N, D<:Number, M} <: AbstractLattice{N, D, M}
    name::Symbol
    coordinates::Matrix{D}
    vectors::SVector{M, SVector{N, D}}
    function Lattice(name::Symbol, coordinates::AbstractMatrix{<:Number}, vectors::SVector{M, <:SVector{N, <:Number}}) where {N, M}
        @assert N==size(coordinates, 1) "Lattice error: shape mismatched."
        datatype = promote_type(Float, eltype(coordinates), eltype(eltype(vectors)))
        coordinates = convert(Matrix{datatype}, coordinates)
        vectors = convert(SVector{M, SVector{N, datatype}}, vectors)
        new{N, datatype, M}(name, coordinates, vectors)
    end
end



"""
    Lattice(coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    Lattice(coordinates::AbstractVector{<:Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing)

Construct a lattice.
"""
function Lattice(coordinates::NTuple{N, Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing) where N
    vectors = isnothing(vectors) ? SVector{0, SVector{N, eltype(eltype(coordinates))}}() : vectorconvert(vectors)
    coordinates = [coordinates[j][i] for i=1:N, j=1:length(coordinates)]
    return Lattice(name, coordinates, vectors)
end

function Lattice(coordinates::AbstractVector{<:Number}...; name::Symbol=:lattice, vectors::Union{AbstractVector{<:AbstractVector{<:Number}}, Nothing}=nothing)
    coordinates = hcat(coordinates...)
    vectors = isnothing(vectors) ? SVector{0, SVector{size(coordinates)[1], eltype(coordinates)}}() : vectorconvert(vectors)
    return Lattice(name, coordinates, vectors)
end

@inline vectorconvert(vectors::SVector{N, <:SVector}) where N = vectors
@inline vectorconvert(vectors::AbstractVector{<:SVector}) = convert(SVector{length(vectors), SVector{length(eltype(vectors)), eltype(eltype(vectors))}}, vectors)
@inline vectorconvert(vectors::AbstractVector{<:AbstractVector}) = convert(SVector{length(vectors), SVector{length(first(vectors)), eltype(eltype(vectors))}}, vectors)



"""
    Point{N, D<:Number}

A point in a unitcell-described lattice.
"""
struct Point{N, D<:Number}
    site::Int
    rcoordinate::SVector{N, D}
    icoordinate::SVector{N, D}
end
@inline Base.:(==)(point₁::Point, point₂::Point) = ==(efficientoperations, point₁, point₂)
@inline Base.isequal(point₁::Point, point₂::Point) = isequal(efficientoperations, point₁, point₂)
@inline Base.show(io::IO, p::Point) = @printf io "Point(%s, %s, %s)" p.site p.rcoordinate p.icoordinate

"""
    Point(site::Integer, rcoordinate::SVector{N, D}, icoordinate::SVector{N, D}) where {N, D<:Number}
    Point(site::Integer, rcoordinate::NTuple{N, <:Number}, icoordinate::NTuple{N, <:Number}=ntuple(i->0, N)) where N
    Point(site::Integer, rcoordinate::AbstractVector{<:Number}, icoordinate::AbstractVector{<:Number}=zero(SVector{length(rcoordinate), Int}))

Construct a labeled point.
"""
@inline function Point(site::Integer, rcoordinate::NTuple{N, <:Number}, icoordinate::NTuple{N, <:Number}=ntuple(i->0, N)) where N
    datatype = promote_type(eltype(rcoordinate), eltype(icoordinate))
    return Point(site, convert(SVector{N, datatype}, rcoordinate), convert(SVector{N, datatype}, icoordinate))
end
@inline function Point(site::Integer, rcoordinate::AbstractVector{<:Number}, icoordinate::AbstractVector{<:Number}=zero(SVector{length(rcoordinate), eltype(rcoordinate)}))
    datatype = promote_type(eltype(rcoordinate), eltype(icoordinate))
    return Point(site, convert(SVector{length(rcoordinate), datatype}, rcoordinate), convert(SVector{length(icoordinate), datatype}, icoordinate))
end

"""
    dimension(point::Point) -> Int
    dimension(::Type{<:Point{N}}) where N -> Int

Get the spatial dimension of a point.
"""
@inline dimension(point::Point) = dimension(typeof(point))
@inline dimension(::Type{<:Point{N}}) where N = N

"""
    dtype(point::Point)
    dtype(::Type{<:Point{N, D} where N}) where {D<:Number}

Get the data type of the coordinates of a point.
"""
@inline dtype(point::Point) = dtype(typeof(point))
@inline dtype(::Type{<:Point{N, D} where N}) where {D<:Number} = D

"""
    isintracell(point::Point) -> Bool

Judge whether a point is intra the unitcell.
"""
@inline isintracell(point::Point) = isapprox(norm(point.icoordinate), 0.0, atol=atol, rtol=rtol)







"""
    Lattice(lattice::AbstractLattice, ranges::NTuple{N, Int}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N)); mode::Symbol=:nonnegative) where N
    Lattice(lattice::AbstractLattice, ranges::NTuple{N, UnitRange{Int}}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N))) where N

Construct a lattice from the translations of another.
"""
function Lattice(lattice::AbstractLattice, ranges::NTuple{N, Int}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N)); mode::Symbol=:nonnegative) where N
    @assert mode∈(:center, :nonnegative) "Lattice error: wrong mode($(repr(mode)))."
    return Lattice(lattice, mode==:center ? map(i->-floor(Int, (i-1)/2):-floor(Int, (i-1)/2)+i-1, ranges) : map(i->0:i-1, ranges), boundaries)
end

function Lattice(lattice::AbstractLattice, ranges::NTuple{N, UnitRange{Int}}, boundaries::NTuple{N, Char}=ntuple(i->'O', Val(N))) where N
    @assert all(map(in(('P', 'O', 'p', 'o')), boundaries)) "Lattice error: boundary conditions must be either 'P'/'p' for 'periodic' or 'O'/'o' for 'open'."
    @assert length(boundaries)==N "Lattice error: mismatched number of ranges and boundaries conditions."
    boundaries = map(uppercase, boundaries)
    name = Symbol(@sprintf "%s%s" getcontent(lattice, :name) join([@sprintf("%s%s%s", boundary=='P' ? "[" : "(", range, boundary=='P' ? "]" : ")") for (range, boundary) in zip(ranges, boundaries)]))
    coordinates = tile(getcontent(lattice, :coordinates), getcontent(lattice, :vectors), product(ranges...))
    vectors = SVector{dimension(lattice), dtype(lattice)}[]
    for (i, vector) in enumerate(getcontent(lattice, :vectors))
        boundaries[i]=='P' && push!(vectors, vector*length(ranges[i]))
    end
    return Lattice(name, coordinates, vectorconvert(vectors))
end





# Generic quantum operator
"""
    QuantumOperator

The abstract type of any quantum operator.
"""
abstract type QuantumOperator end
@inline Base.:(==)(m₁::QuantumOperator, m₂::QuantumOperator) = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::QuantumOperator, m₂::QuantumOperator) = isequal(efficientoperations, m₁, m₂)







# Operator unit
"""
    OperatorUnit <: QuantumOperator

An operator unit is the irreducible symbolic unit to represent a quantum operator.

It plays the role of the symbols as in usual computer algebras while it can host internal structures, which is convenient for quantum operators in representative of the internal degrees of freedom.
"""
abstract type OperatorUnit <: QuantumOperator end
@inline Base.show(io::IO, u::OperatorUnit) = @printf io "%s(%s)" nameof(typeof(u)) join(map(repr, ntuple(i->getfield(u, i), Val(fieldcount(typeof(u))))), ", ")
@inline @generated Base.hash(u::OperatorUnit, h::UInt) = Expr(:call, :hash, Expr(:tuple, [:(getfield(u, $i)) for i=1:fieldcount(u)]...), :h)






# ID of a composite quantum operator
"""
    ID{U<:OperatorUnit, N}

The id of a composite quantum operator, which is an ordered set of operator units.

Type alias for `NTuple{N, U} where {U<:OperatorUnit}`.
"""
const ID{U<:OperatorUnit, N} = NTuple{N, U}
@inline Base.promote_rule(::Type{Tuple{}}, I::Type{<:Tuple{OperatorUnit, Vararg{OperatorUnit}}}) = ID{I|>eltype}
@inline Base.promote_rule(I::Type{<:Tuple{OperatorUnit, Vararg{OperatorUnit}}}, ::Type{Tuple{}}) = ID{I|>eltype}

"""
    ID(id::OperatorUnit...)
    ID(u::OperatorUnit, id::ID{OperatorUnit})
    ID(id::ID{OperatorUnit}, u::OperatorUnit)
    ID(id₁::ID{OperatorUnit}, id₂::ID{OperatorUnit})

Get the id from operator units/ids.
"""
@inline ID(id::OperatorUnit...) = id
@inline ID(u::OperatorUnit, id::ID{OperatorUnit}) = ID(u, id...)
@inline ID(id::ID{OperatorUnit}, u::OperatorUnit) = ID(id..., u)
@inline ID(id₁::ID{OperatorUnit}, id₂::ID{OperatorUnit}) = ID(id₁..., id₂...)

"""
    ID(::Type{U}, attrs::Vararg{NTuple{N}, M}) where {U<:OperatorUnit, N, M}

Get the composite id from the components of singular ids.
"""
@inline @generated function ID(::Type{U}, attrs::Vararg{NTuple{N, Any}, M}) where {U<:OperatorUnit, N, M}
    exprs = []
    for i = 1:N
        args = [:(attrs[$j][$i]) for j = 1:M]
        push!(exprs, :(U($(args...))))
    end
    return :(ID($(exprs...)))
end









# IID and Internal
"""
    IID <: OperatorUnit

The id of an internal degree of freedom.
"""
abstract type IID <: OperatorUnit end






# Vector spaces
"""
    VectorSpace{B} <: AbstractVector{B}

Abstract vector space.
"""
abstract type VectorSpace{B} <: AbstractVector{B} end
@inline Base.:(==)(vs₁::VectorSpace, vs₂::VectorSpace) = ==(efficientoperations, vs₁, vs₂)
@inline Base.isequal(vs₁::VectorSpace, vs₂::VectorSpace) = isequal(efficientoperations, vs₁, vs₂)
@inline Base.size(vs::VectorSpace) = (length(vs),)








"""
    Internal{I<:IID} <: VectorSpace{I}

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: VectorSpace{I} end







"""
    CompositeDict{K, V}

A composite dict can be considered as a dict that is implemented by including a concrete subtype of `AbstractDict` as its data attribute.
"""
abstract type CompositeDict{K, V} <: AbstractDict{K, V} end
@inline contentnames(::Type{<:CompositeDict}) = (:contents,)
@inline dissolve(cd::CompositeDict, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(cd, :contents), args...; kwargs...)

@inline Base.isempty(cd::CompositeDict) = isempty(getcontent(cd, :contents))
@inline Base.length(cd::CompositeDict) = length(getcontent(cd, :contents))
@inline Base.haskey(cd::CompositeDict, key) = haskey(getcontent(cd, :contents), key)
@inline Base.in(p::Pair, cd::CompositeDict, valcmp=(==)) = in(p, getcontent(cd, :contents), valcmp)
@inline Base.:(==)(cd1::CompositeDict, cd2::CompositeDict) = ==(efficientoperations, cd1, cd2)
@inline Base.isequal(cd1::CompositeDict, cd2::CompositeDict) = isequal(efficientoperations, cd1, cd2)
@inline Base.get(cd::CompositeDict, key, default) = get(getcontent(cd, :contents), key, default)
@inline Base.get(f::Base.Callable, cd::CompositeDict, key) = get(f, getcontent(cd, :contents), key)
@inline Base.getkey(cd::CompositeDict, key, default) = getkey(getcontent(cd, :contents), key, default)
@inline Base.getindex(cd::CompositeDict{K, V}, index::K) where {K, V} = getcontent(cd, :contents)[index]
@inline Base.push!(cd::CompositeDict, ps::Pair...) = (push!(getcontent(cd, :contents), ps...); cd)
@inline Base.get!(cd::CompositeDict, key, default) = get!(getcontent(cd, :contents), key, default)
@inline Base.get!(default::Union{Function, Type}, cd::CompositeDict, key) = get!(default, getcontent(cd, :contents), key)
@inline Base.setindex!(cd::CompositeDict{K, V}, value::V, index::K) where {K, V} = (getcontent(cd, :contents)[index] = value)
@inline Base.pop!(cd::CompositeDict) = pop!(getcontent(cd, :contents))
@inline Base.pop!(cd::CompositeDict, key) = pop!(getcontent(cd, :contents), key)
@inline Base.pop!(cd::CompositeDict, key, default) = pop!(getcontent(cd, :contents), key, default)
@inline Base.delete!(cd::CompositeDict, key) = (delete!(getcontent(cd, :contents), key); cd)
@inline Base.empty!(cd::CompositeDict) = (empty!(getcontent(cd, :contents)); cd)
@inline Base.merge(cd::CD, others::CD...) where CD <: CompositeDict = merge!(empty(cd), cd, others...)
@inline Base.merge(combine::Function, cd::CD, others::CD...) where {CD<:CompositeDict} = merge!(combine, empty(cd), cd, others...)
@inline Base.empty(cd::CompositeDict) = rawtype(typeof(cd))(dissolve(cd, empty)...)
@inline Base.iterate(cd::CompositeDict) = iterate(getcontent(cd, :contents))
@inline Base.iterate(cd::CompositeDict, state) = iterate(getcontent(cd, :contents), state)
@inline Base.keys(cd::CompositeDict) = keys(getcontent(cd, :contents))
@inline Base.values(cd::CompositeDict) = values(getcontent(cd, :contents))
@inline Base.pairs(cd::CompositeDict) = pairs(getcontent(cd, :contents))






# Index and CompositeIndex
"""
    Index{S<:Union{Int, Colon}, I<:SimpleIID} <: OperatorUnit

The index of a degree of freedom, which consist of the spatial part and the internal part.
"""
struct Index{S<:Union{Int, Colon}, I<:SimpleIID} <: OperatorUnit
    site::S
    iid::I
end
@inline parameternames(::Type{<:Index}) = (:site, :iid)
@inline isparameterbound(::Type{<:Index}, ::Val{:site}, ::Type{S}) where {S<:Union{Int, Colon}} = !isconcretetype(S)
@inline isparameterbound(::Type{<:Index}, ::Val{:iid}, ::Type{I}) where {I<:SimpleIID} = !isconcretetype(I)
@inline isdefinite(index::Index) = isdefinite(typeof(index))
@inline isdefinite(::Type{<:Index{<:Union{Int, Colon}, I}}) where {I<:SimpleIID} = isdefinite(I)
@inline isdefinite(indexes::Tuple{Vararg{Index}}) = isdefinite(typeof(indexes))
@inline @generated isdefinite(::Type{T}) where {T<:Tuple{Vararg{Index}}} = Expr(:call, :all, Expr(:tuple, [:(isdefinite(fieldtype(T, $i))) for i=1:fieldcount(T)]...))
@inline Base.show(io::IO, index::Index{Colon}) = @printf io "Index(:, %s)" index.iid

"""
    iidtype(index::Index)
    iidtype(::Type{I}) where {I<:Index}

Get the type of the internal part of an index.
"""
@inline iidtype(index::Index) = iidtype(typeof(index))
@inline iidtype(::Type{I}) where {I<:Index} = parametertype(I, :iid)

"""
    statistics(index::Index) -> Symbol
    statistics(::Type{<:Index{I}}) where {I<:SimpleIID} -> Symbol

Get the statistics of an index.
"""
@inline statistics(index::Index) = statistics(typeof(index))
@inline statistics(::Type{<:Index{<:Union{Int, Colon}, I}}) where {I<:SimpleIID} = statistics(I)




"""
    AbstractCompositeIndex{I<:Index} <: OperatorUnit

The abstract type of a composite index.
"""
abstract type AbstractCompositeIndex{I<:Index} <: OperatorUnit end
@inline contentnames(::Type{<:AbstractCompositeIndex}) = (:index,)
@inline parameternames(::Type{<:AbstractCompositeIndex}) = (:index,)
@inline isparameterbound(::Type{<:AbstractCompositeIndex}, ::Val{:index}, ::Type{I}) where {I<:Index} = !isconcretetype(I)

"""
    indextype(::AbstractCompositeIndex)
    indextype(::Type{<:AbstractCompositeIndex})

Get the index type of a composite index.
"""
@inline indextype(index::AbstractCompositeIndex) = indextype(typeof(index))
@inline @generated indextype(::Type{I}) where {I<:AbstractCompositeIndex} = parametertype(supertype(I, :AbstractCompositeIndex), :index)

"""
    statistics(index::AbstractCompositeIndex) -> Symbol
    statistics(::Type{<:AbstractCompositeIndex{I}}) where {I<:Index} -> Symbol

Get the statistics of a composite operator id.
"""
@inline statistics(index::AbstractCompositeIndex) = statistics(typeof(index))
@inline statistics(::Type{<:AbstractCompositeIndex{I}}) where {I<:Index} = statistics(I)

"""
    CompositeIndex{I<:Index, V<:SVector} <: AbstractCompositeIndex{I}

Composite index of a quantum operator.
"""
struct CompositeIndex{I<:Index, V<:SVector} <: AbstractCompositeIndex{I}
    index::I
    rcoordinate::V
    icoordinate::V
    CompositeIndex(index::Index, rcoordinate::V, icoordinate::V) where {V<:SVector} = new{typeof(index), V}(index, compositeindexcoordinate(rcoordinate), compositeindexcoordinate(icoordinate))
end
@inline contentnames(::Type{<:CompositeIndex}) = (:index, :rcoordinate, :icoordinate)
@inline parameternames(::Type{<:CompositeIndex}) = (:index, :coordination)
@inline isparameterbound(::Type{<:CompositeIndex}, ::Val{:coordination}, ::Type{V}) where {V<:SVector} = !isconcretetype(V)
@inline Base.hash(index::CompositeIndex, h::UInt) = hash((index.index, Tuple(index.rcoordinate)), h)
@inline Base.propertynames(::ID{CompositeIndex}) = (:indexes, :rcoordinates, :icoordinates)
@inline Base.show(io::IO, index::CompositeIndex) = @printf io "CompositeIndex(%s, %s, %s)" index.index index.rcoordinate index.icoordinate
@inline compositeindexcoordinate(vector::SVector) = vector
@inline compositeindexcoordinate(vector::SVector{N, Float}) where N = SVector(ntuple(i->vector[i]===-0.0 ? 0.0 : vector[i], Val(N)))

"""
    CompositeIndex(index::Index, rcoordinate, icoordinate)
    CompositeIndex(index::Index; rcoordinate, icoordinate)

Construct an operator id.
"""
@inline CompositeIndex(index::Index, rcoordinate, icoordinate) = CompositeIndex(index, SVector{length(rcoordinate)}(rcoordinate), SVector{length(icoordinate)}(icoordinate))
@inline CompositeIndex(index::Index; rcoordinate, icoordinate) = CompositeIndex(index, rcoordinate, icoordinate)

"""
    adjoint(index::CompositeIndex) -> typeof(index)

Get the adjoint of an operator id.
"""
@inline Base.adjoint(index::CompositeIndex) = CompositeIndex(index.index', index.rcoordinate, index.icoordinate)

"""
    indextype(I::Type{<:SimpleInternal}, P::Type{<:Point}, ::Val)

Get the compatible composite index type based on the information of its internal part.
"""
@inline function indextype(I::Type{<:SimpleInternal}, P::Type{<:Point}, ::Val)
    return fulltype(CompositeIndex, NamedTuple{(:index, :coordination), Tuple{fulltype(Index, NamedTuple{(:site, :iid), Tuple{Int, eltype(I)}}), SVector{dimension(P), dtype(P)}}})
end

"""
    rcoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}}) -> SVector

Get the whole rcoordinate of an operator.
"""
@inline function rcoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}})
    rank(opt)==1 && return id(opt)[1].rcoordinate
    rank(opt)==2 && return id(opt)[2].rcoordinate-id(opt)[1].rcoordinate
    error("rcoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}}) -> SVector

Get the whole icoordinate of an operator.
"""
@inline function icoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}})
    rank(opt)==1 && return id(opt)[1].icoordinate
    rank(opt)==2 && return id(opt)[2].icoordinate-id(opt)[1].icoordinate
    error("icoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    script(::Val{:rcoordinate}, index::CompositeIndex; kwargs...) -> String
    script(::Val{:icoordinate}, index::CompositeIndex; kwargs...) -> String

Get the `rcoordinate/icoordinate` script of a composite index.
"""
@inline script(::Val{:rcoordinate}, index::CompositeIndex; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.rcoordinate), ", ")
@inline script(::Val{:icoordinate}, index::CompositeIndex; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.icoordinate), ", ")

"""
    script(::Val{:integercoordinate}, index::CompositeIndex; vectors, kwargs...)

Get the integral script of the icoordinate of an composite index.
"""
function script(::Val{:integercoordinate}, index::CompositeIndex; vectors, kwargs...)
    rcoeff = decompose(index.icoordinate, vectors...)
    icoeff = Int.(round.(rcoeff))
    @assert isapprox(efficientoperations, rcoeff, icoeff) "script error: mismatched icoordinate of the input composite index and vectors."
    return @sprintf "[%s]" join(icoeff, ", ")
end

"""
    script(::Val{attr}, index::CompositeIndex; kwargs...) where attr

Get the `attr` script of an index, which is contained in its index.
"""
@inline script(::Val{attr}, index::CompositeIndex; kwargs...) where attr = script(Val(attr), index.index; kwargs...)

"""
    script(::Val{:site}, index::Index; kwargs...) -> String
    script(attr::Val, index::Index; kwargs...) -> String

Get the required script of a spin index.
"""
@inline script(::Val{:site}, index::Index; kwargs...) = string(index.site)
@inline script(::Val{:site}, index::Index{Colon}; kwargs...) = ":"
@inline script(attr::Val, index::Index; kwargs...) = script(attr, index.iid; kwargs...)








# Hilbert
"""
    Hilbert{I<:Internal} <: CompositeDict{Int, I}

Hilbert space at a lattice.
"""
struct Hilbert{I<:Internal} <: CompositeDict{Int, I}
    contents::Dict{Int, I}
    Hilbert(contents::Dict{Int, <:Internal}) = new{valtype(contents)}(contents)
end

"""
    Hilbert(ps::Pair...)
    Hilbert(kv)

Construct a Hilbert space the same way as a Dict.
"""
@inline Hilbert(ps::Pair...) = Hilbert(ps)
@inline Hilbert(kv) = Hilbert(Dict(kv))

"""
    Hilbert(internals::Internal...)
    Hilbert(internals::Tuple{Vararg{Internal}})
    Hilbert(internals::AbstractVector{<:Internal})

Construct a Hilbert space with the given internals.
"""
@inline Hilbert(internals::Internal...) = Hilbert(internals)
@inline Hilbert(internals::Tuple{Vararg{Internal}}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))
@inline Hilbert(internals::AbstractVector{<:Internal}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))

"""
    Hilbert(internal::Internal, num::Int)

Construct a Hilbert space with all internal spaces the same.
"""
@inline Hilbert(internal::Internal, num::Int) = Hilbert(i=>internal for i in 1:num)







"""
    SimpleIID <: IID

The id of a simple internal degree of freedom.
"""
abstract type SimpleIID <: IID end
@inline statistics(iid::SimpleIID) = statistics(typeof(iid))
@inline isdefinite(iid::SimpleIID) = isdefinite(typeof(iid))
@inline isdefinite(::Type{<:SimpleIID}) = false







"""
    SimpleInternal{I<:SimpleIID} <: Internal{I}

The simple internal degrees of freedom at a single point.
"""
abstract type SimpleInternal{I<:SimpleIID} <: Internal{I} end
@inline VectorSpaceStyle(::Type{<:SimpleInternal}) = VectorSpaceCartesian()
@inline Base.show(io::IO, i::SimpleInternal) = @printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i, name))" for name in i|>typeof|>fieldnames), ", ")






# Canonical complex fermionic/bosonic systems and hardcore bosonic systems
## FID
"""
    annihilation

Indicate that the nambu index is annihilation.
"""
const annihilation = 1

"""
    creation

Indicate that the nambu index is creation.
"""
const creation = 2

@inline default(::Colon) = ":"
@inline default(value::Char) = repr(value)
@inline default(value) = string(value)
@inline default(value::Rational{Int}) = value.den==1 ? repr(value.num) : repr(value)






"""
    FID{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} <: SimpleIID

The Fock id.
"""
struct FID{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} <: SimpleIID
    orbital::O
    spin::S
    nambu::N
    function FID{T}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) where T
        @assert T∈(:f, :b, wildcard) "FID error: wrong statistics."
        isa(spin, Rational{Int}) && @assert spin.den∈(1, 2) "FID error: wrong spin."
        isa(nambu, Int) && @assert nambu∈(1, 2) "FID error: wrong input nambu($nambu)."
        new{T, typeof(orbital), typeof(spin), typeof(nambu)}(orbital, spin, nambu)
    end
end
@inline FID{T}(orbital::Union{Int, Symbol, Colon}, spin::Int, nambu::Union{Int, Symbol, Colon}) where T = FID{T}(orbital, spin//1, nambu)
@inline FID(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) = FID{wildcard}(orbital, spin, nambu)
@inline Base.:(==)(fid₁::FID, fid₂::FID) = statistics(fid₁)==statistics(fid₂) && ==(efficientoperations, fid₁, fid₂)
@inline Base.isequal(fid₁::FID, fid₂::FID) = isequal(statistics(fid₁), statistics(fid₂)) && isequal(efficientoperations, fid₁, fid₂)
@inline Base.hash(fid::FID, h::UInt) = hash((statistics(fid), fid.orbital, fid.spin, fid.nambu), h)
@inline Base.show(io::IO, fid::FID) = @printf io "FID{%s}(%s)" repr(statistics(fid)) join((fid.orbital|>default, fid.spin|>default, fid.nambu|>default), ", ")
@inline Base.show(io::IO, fid::FID{wildcard}) = @printf io "FID(%s)" join((fid.orbital|>default, fid.spin|>default, fid.nambu|>default), ", ")
@inline Base.adjoint(fid::FID{T, <:Union{Int, Symbol, Colon}, <:Union{Rational{Int}, Symbol, Colon}, Int}) where T = FID{T}(fid.orbital, fid.spin, 3-fid.nambu)
@inline @generated function Base.replace(fid::FID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(fid, $name))) for name in QuoteNode.(fieldnames(fid))]
    return :(rawtype(typeof(fid)){statistics(fid)}($(exprs...)))
end
@inline statistics(::Type{<:FID{T}}) where T = T

#@inline FID(iid::FID, ::CompositeInternal) = iid







## Fock
"""
    Fock{T} <: SimpleInternal{FID{T, Int, Rational{Int}, Int}}

The Fock internal degrees of freedom.
"""
struct Fock{T} <: SimpleInternal{FID{T, Int, Rational{Int}, Int}}
    norbital::Int
    nspin::Int
    function Fock{T}(norbital::Int, nspin::Int) where T
        @assert T∈(:f, :b) "Fock error: wrong statistics."
        new{T}(norbital, nspin)
    end
end
@inline Base.eltype(::Type{Fock}) = (FID{T, Int, Rational{Int}, Int} where T)
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, 1:2)
@inline Base.CartesianIndex(fid::FID{T}, fock::Fock{T}) where T = CartesianIndex(fid.orbital, Int(fid.spin+(fock.nspin-1)//2)+1, fid.nambu)

@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2]-1-(fock.nspin-1)//2, index[3])

@inline Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))

@inline Base.show(io::IO, fock::Fock) = @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")

@inline Base.match(::Type{<:FID{wildcard}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T₁}}, ::Type{<:Fock{T₂}}) where {T₁, T₂} = false









# Combinatorics
"""
    Combinatorics{M, C}

Abstract combinatorial algorithms.
"""
abstract type Combinatorics{M, C} end
@inline Base.eltype(::Type{<:Combinatorics{M, C}}) where {M, C} = NTuple{M, eltype(C)}




"""
    Combinations{M}(contents::C) where {M, C}

Combinations of M elements from contents. Duplicates are not allowed.
"""
struct Combinations{M, C} <: Combinatorics{M, C}
    contents::C
    N::Int
    Combinations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(c::Combinations{M}) where M = binomial(c.N, M)
Base.iterate(c::Combinations{M}) where M = (M > c.N) ? nothing : (M == 0) ? ((), [c.N+2]) : (ntuple(i->c.contents[i], Val(M)), nextmstate!(collect(1:M), c.N, M))
Base.iterate(c::Combinations{M}, state) where M = (state[1] > c.N-M+1) ? nothing : (ntuple(i->c.contents[state[i]], Val(M)), nextmstate!(state, c.N, M))
function nextmstate!(state::Vector{Int}, N::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        (state[i] > N-(M-i)) && continue
        for j = i+1:M
            state[j] = state[j-1] + 1
        end
        break
    end
    state
end













"""
    Diagonal{Fields} <: Function

Construct a pattern for a set of homogenous `Index`es that all the specified fields of their contained iids should be diagonal, respectively.
"""
struct Diagonal{Fields} <: Function
    Diagonal(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline Diagonal(fields::Symbol...) = Diagonal(fields)
@inline Diagonal(::Type{I}) where {I<:Index} = _diagonal(I, I|>diagonalizablefields|>Val)
@generated function _diagonal(::Type{I}, ::Val{choices}) where {I<:Index, choices}
    fields = Symbol[]
    T = iidtype(I)
    for choice in choices
        fieldtype(T, choice)==Colon && push!(fields, choice)
    end
    return Diagonal(fields...)
end
@inline @generated diagonalizablefields(::Type{I}) where {I<:Index} = fieldnames(iidtype(I))
@generated function (diagonal::Diagonal{fields})(indexes::Tuple{I, Vararg{I, N}}) where {fields, I<:Index, N}
    exprs = []
    for field in QuoteNode.(fields)
        push!(exprs, Expr(:call, allequal, Expr(:tuple, [:(getfield(indexes[$i].iid, $field)) for i=1:(N+1)]...)))
    end
    return Expr(:call, all, Expr(:tuple, exprs...))
end







"""
    Constraint{RS, N, C<:NTuple{N, Function}}

The constraint of the indexes of internal degrees of freedom in a coupling.
"""
struct Constraint{RS, N, C<:NTuple{N, Function}}
    representations::NTuple{N, String}
    conditions::C
    function Constraint{RS}(representations::NTuple{N, String}, conditions::NTuple{N, Function})  where {RS, N}
        @assert isa(RS, NTuple) "Constraint error: ranks (`RS`) must be tuple of integers."
        @assert length(RS)==N "Constraint error: mismatched number of ranks ($RS), representations ($representations) and conditions ($conditions)."
        new{RS, N, typeof(conditions)}(representations, conditions)
    end
end
@inline Base.:(==)(constraint₁::Constraint{RS₁}, constraint₂::Constraint{RS₂}) where {RS₁, RS₂} = RS₁==RS₂ && constraint₁.representations==constraint₂.representations
@inline Base.:isequal(constraint₁::Constraint{RS₁}, constraint₂::Constraint{RS₂}) where {RS₁, RS₂} = isequal(RS₁, RS₂) && isequal(constraint₁.representations, constraint₂.representations)
@inline Base.hash(constraint::Constraint{RS}, h::UInt) where RS = hash((RS, constraint.representations), h)


"""
    Constraint{R}() where R
    Constraint{R}(condition::Union{Pattern, Diagonal}) where R
    Constraint{R}(representation::String, condition::Function) where R

Construct a constraint with only one condition.
"""
@inline Constraint{R}() where R = Constraint{R}(Diagonal())
#@inline Constraint{R}(condition::Union{Pattern, Diagonal}) where R = Constraint{R}("pattern", condition)
# ffr
@inline Constraint{R}(representation::String, condition::Function) where R = Constraint{(R,)}((representation,), (condition,))

"""
    Constraint(indexes::Index...)
    Constraint(indexes::NTuple{N, Index}) where N

Construct a constraint based on the pattern of the input indexes.
"""
@inline Constraint(index::Index, indexes::Index...) = Constraint((index, indexes...))
@inline Constraint(indexes::Tuple{Index, Vararg{Index}}) = Constraint{fieldcount(typeof(indexes))}("pattern", constraint(indexes, indexes|>typeof|>isdiagonalable|>Val))
@inline constraint(indexes::Tuple{Index, Vararg{Index}}, ::Val{true}) = Diagonal(eltype(indexes))
@inline constraint(indexes::Tuple{Index, Vararg{Index}}, ::Val{false}) = Pattern(indexes)
@generated function isdiagonalable(::Type{T}) where {T<:Tuple{Index, Vararg{Index}}}
    types = fieldtypes(T)
    allequal(types) || return false
    for type in map(iidtype, types)
        for i = 1:fieldcount(type)
            fieldtype(type, i)==Symbol && return false
        end
    end
    return true
end

"""
    rank(constraint::Constraint) -> Int
    rank(::Type{<:Constraint{RS}}) where RS -> Int

Get the rank of the coupling indexes that a constraint can apply.
"""
@inline rank(constraint::Constraint) = rank(typeof(constraint))
@inline @generated rank(::Type{<:Constraint{RS}}) where RS = sum(RS)

"""
    rank(constraint::Constraint, i::Integer) -> Int
    rank(::Type{<:Constraint{RS}}, i::Integer) where RS -> Int

Get the rank of the ith homogenous segment of the coupling indexes that a constraint can apply.
"""
@inline rank(constraint::Constraint, i::Integer) = rank(typeof(constraint), i)
@inline rank(::Type{<:Constraint{RS}}, i::Integer) where RS = RS[i]

"""
    match(constraint::Constraint, indexes::Tuple{Vararg{Index}}) -> Bool

Judge whether a composite iid fulfills a constraint.
"""
@generated function Base.match(constraint::Constraint{RS}, indexes::Tuple{Index, Vararg{Index}}) where RS
    @assert rank(constraint)==fieldcount(indexes) "match error: mismatched rank of indexes and constraint."
    exprs, count = [], 1
    for (i, r) in enumerate(RS)
        start, stop = count, count+r-1
        segment = Expr(:tuple, [:(indexes[$pos]) for pos=start:stop]...)
        push!(exprs, :(constraint.conditions[$i]($segment)::Bool))
        count = stop+1
    end
    return Expr(:call, :all, Expr(:tuple, exprs...))
end

"""
    *(constraint₁::Constraint, constraint₂::Constraint) -> Constraint

Get the combination of two sets of constraints.
"""
@inline function Base.:*(constraint₁::Constraint{RS₁}, constraint₂::Constraint{RS₂}) where {RS₁, RS₂}
    return Constraint{(RS₁..., RS₂...)}((constraint₁.representations..., constraint₂.representations...), (constraint₁.conditions..., constraint₂.conditions...))
end

"""
    @indexes index₁ index₂ ...
    @indexes(index₁, index₂, ...; constraint=...)

Construct an set of indexes and its constraint according to the input index pattern and an optional constraint.
"""
macro indexes(exprs...)
    if exprs[1].head==:parameters
        @assert(
            length(exprs[1].args)==1 && exprs[1].args[1].head==:kw && exprs[1].args[1].args[1]==:constraint,
            "@indexes error: constraint must be specified by the keyword argument `constraint`."
        )
        constraint = exprs[1].args[1].args[2]
        patterns = exprs[2:end]
    else
        constraint = true
        patterns = exprs
    end
    indexes = []
    blocks = []
    for (i, pattern) in enumerate(patterns)
        @assert pattern.head==:call && pattern.args[end].head==:call "@indexes error: wrong pattern."
        attrs = []
        for (j, attr) in enumerate(pattern.args[end].args[2:end])
            isa(attr, Expr) && !(attr.head==:call && length(attr.args)==3 && attr.args[1]==:// && isa(attr.args[2], Int) && isa(attr.args[3], Int)) && (attr = Symbol(attr))
            if isa(attr, Symbol) && attr≠wildcard
                @assert Base.isidentifier(attr) "@indexes error: wrong pattern."
                push!(blocks, quote
                    local $attr
                    if @isdefined($attr)
                        ($attr)!=getfield(indexes[$i].iid, $j) && return false
                    else
                        $attr = getfield(indexes[$i].iid, $j)
                    end
                end)
                attr = QuoteNode(attr)
            end
            push!(attrs, attr)
        end
        push!(indexes, Expr(:call, pattern.args[1], pattern.args[2], Expr(:call, pattern.args[end].args[1], attrs...)))
    end
    N = length(indexes)
    indexes = Expr(:tuple, indexes...)
    blocks = Expr(:block, vcat([block.args for block in blocks]...)...)
    name = gensym("constraint")
    representation = constraint==true ? "pattern" : string(constraint)
    return quote
        function ($name)(indexes::NTuple{$N, Index})
            $blocks
            return $constraint
        end
        ($(esc(indexes)), Constraint{$N}($representation, $name))
    end
end











# Operator pack
"""
    OperatorPack{V, I<:Tuple} <: QuantumOperator

The entity that represent the pack of a number and several quantum units.

Basically, a concrete subtype should contain two predefined contents:
- `value::V`: the coefficient of the pack
- `id::I`: the total id of the pack
"""
abstract type OperatorPack{V, I<:Tuple} <: QuantumOperator end
@inline contentnames(::Type{<:OperatorPack}) = (:value, :id)
@inline parameternames(::Type{<:OperatorPack}) = (:value, :id)
@inline isparameterbound(::Type{<:OperatorPack}, ::Val{:value}, ::Type{V}) where V = false
@inline isparameterbound(::Type{<:OperatorPack}, ::Val{:id}, ::Type{I}) where {I<:Tuple} = !isconcretetype(I)
@inline function Base.promote_rule(::Type{M₁}, ::Type{M₂}) where {M₁<:OperatorPack, M₂<:OperatorPack}
    M₁<:M₂ && return M₂
    M₂<:M₁ && return M₁
    r₁, r₂ = rank(M₁), rank(M₂)
    M = r₂==0 ? rawtype(M₁) : r₁==0 ? rawtype(M₂) : typejoin(rawtype(M₁), rawtype(M₂))
    return fulltype(M, promoteparameters(parameterpairs(M₁), parameterpairs(M₂)))
end
@inline Base.promote_rule(::Type{M}, ::Type{N}) where {M<:OperatorPack, N<:Number} = reparameter(M, :value, promote_type(valtype(M), N))

"""
    rank(m::OperatorPack) -> Int
    rank(::Type{M}) where {M<:OperatorPack} -> Int

Get the rank of an `OperatorPack`.
"""
@inline rank(m::OperatorPack) = rank(typeof(m))
@inline rank(::Type{M}) where {M<:OperatorPack} = rank(idtype(M))














## Coupling
"""
    Coupling{V, I<:ID{Index}, C<:Constraint} <: OperatorPack{V, Tuple{I, C}}

The coupling intra/inter internal degrees of freedom at different lattice points.
"""
struct Coupling{V, I<:ID{Index}, C<:Constraint} <: OperatorPack{V, Tuple{I, C}}
    value::V
    indexes::I
    constraint::C
end
@inline parameternames(::Type{<:Coupling}) = (:value, :indexes, :constraint)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:indexes}, ::Type{I}) where {I<:ID{Index}} = !isconcretetype(I)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:constraint}, ::Type{C}) where {C<:Constraint} = !isconcretetype(C)
@inline idtype(M::Type{<:Coupling}) = Tuple{parametertype(M, :indexes), parametertype(M, :constraint)}
@inline getcontent(coupling::Coupling, ::Val{:id}) = (coupling.indexes, coupling.constraint)
@inline rank(M::Type{<:Coupling}) = fieldcount(parametertype(M, :indexes))
@inline Coupling(id::Tuple{ID{Index}, Constraint}) = Coupling(1, id)
@inline Coupling(value, id::Tuple{ID{Index}, Constraint}) = Coupling(value, id...)
@inline CompositeIID(coupling::Coupling, ::Val=Val(:term)) = CompositeIID(coupling.indexes.iids)
@inline Constraint(coupling::Coupling, ::Val=Val(:term)) = coupling.constraint
@inline Base.iterate(coupling::Coupling) = (coupling, nothing)
@inline Base.iterate(coupling::Coupling, ::Nothing) = nothing
@inline Base.eltype(coupling::Coupling) = eltype(typeof(coupling))
@inline Base.eltype(C::Type{<:Coupling}) = C
@inline Base.length(coupling::Coupling) = length(typeof(coupling))
@inline Base.length(::Type{<:Coupling}) = 1
@inline function Base.show(io::IO, coupling::Coupling)
    @printf io "%s" (coupling.value≈1 ? "" : coupling.value≈-1 ? "- " : string(decimaltostr(coupling.value), " "))
    len, count = length(coupling.constraint.representations), 1
    for i = 1:len
        r = rank(coupling.constraint, i)
        start, stop = count, count+r-1
        indexes = coupling.indexes[start:stop]
        if i==1
            @printf io "%s" (isdefinite(indexes) ? "" : "∑")
            (len>1 || !isdefinite(indexes)) && @printf io "%s" "["
        else
            @printf io " ⋅ %s[" (isdefinite(indexes) ? "" : "∑")
        end
        @printf io "%s" join(indexes, " ")
        (len>1 || !isdefinite(indexes)) && @printf io "%s" "]"
        coupling.constraint.representations[i]=="pattern" || @printf io "(%s)" coupling.constraint.representations[i]
        count = stop+1
    end
end
@inline Base.summary(io::IO, couplings::Vector{<:Coupling}) = @printf io "%s-element Vector{Coupling}" length(couplings)




"""
    Coupling(indexes::Index...)
    Coupling(value, indexes::Index...)
    Coupling(value, indexes::Tuple{Vararg{Index}})

Construct a `Coupling` with the input `indexes` as the pattern.
"""
@inline Coupling(index::Index, indexes::Index...) = Coupling(1, (index, indexes...))
@inline Coupling(value, index::Index, indexes::Index...) = Coupling(value, (index, indexes...))
@inline Coupling(value, indexes::Tuple{Index, Vararg{Index}}) = Coupling(value, indexes, Constraint(indexes))

"""
    Coupling(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    Coupling(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    Coupling{N}(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    Coupling{N}(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}

Construct a `Coupling` with the input sites and the fields of a kind of simple iid.
"""


@inline Coupling(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID} = Coupling{N}(sites, I, fields...)
@inline Coupling(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID} = Coupling{N}(value, sites, I, fields...)
@inline Coupling{N}(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID} = Coupling{N}(1, sites, I, fields...)
@inline function Coupling{N}(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    return Coupling(value, ID(Index, default(sites, Val(N)), ID(I, map(field->default(field, Val(N)), fields)...)))
end
@inline default(fields, ::Val) = fields
@inline default(::Colon, N::Val) = ntuple(i->:, N)



"""
    *(cp₁::Coupling, cp₂::Coupling) -> Coupling

Get the multiplication between two coupling.
"""
@inline Base.:*(cp₁::Coupling, cp₂::Coupling) = Coupling(cp₁.value*cp₂.value, (cp₁.indexes..., cp₂.indexes...), cp₁.constraint*cp₂.constraint)





# Operator prod
"""
    OperatorProd{V, I<:ID{OperatorUnit}} <: OperatorPack{V, I}

A special kind of `OperatorPack`, where the relation between the coefficient and quantum units could be viewed as product.
"""
abstract type OperatorProd{V, I<:ID{OperatorUnit}} <: OperatorPack{V, I} end
@inline Base.eltype(m::OperatorProd) = eltype(typeof(m))
@inline Base.eltype(::Type{M}) where {M<:OperatorProd} = eltype(idtype(M))
@inline Base.iterate(m::OperatorProd) = iterate(id(m))
@inline Base.iterate(m::OperatorProd, state) = iterate(id(m), state)

"""
    length(m::OperatorProd) -> Int

Get the length of an `OperatorProd`.
"""
@inline Base.length(m::OperatorProd) = rank(m)
@inline Base.firstindex(m::OperatorProd) = 1
@inline Base.lastindex(m::OperatorProd) = rank(m)

"""
    getindex(m::OperatorProd, i::Integer) -> OperatorUnit
    getindex(m::OperatorProd, slice) -> OperatorProd

Overloaded `[]`.
"""
@inline Base.getindex(m::OperatorProd, i::Integer) = id(m)[i]
@inline Base.getindex(m::OperatorProd, slice) = rawtype(typeof(m))(dissolve(m, getindex, (slice,))...)
@inline dissolve(m::OperatorProd, ::Val{:value}, ::typeof(getindex), ::Tuple{Any}, ::NamedTuple) = one(valtype(m))
@inline dissolve(m::OperatorProd, ::Val{:id}, ::typeof(getindex), slice::Tuple{Any}, ::NamedTuple) = id(m)[first(slice)]

"""
    split(m::OperatorProd) -> Tuple{valtype(m), Vararg{OperatorUnit}}

Split an `OperatorProd` into the coefficient and a sequence of `OperatorUnit`s.
"""
@inline Base.split(m::OperatorProd) = (value(m), id(m)...)

"""
    sequence(m::OperatorProd, table) -> NTuple{rank(m), Int}

Get the sequence of the id of a quantum operator according to a table.
"""
@inline sequence(m::OperatorProd, table) = map(u->table[u], id(m))





# Operator
"""
    Operator{V<:Number, I<:ID{OperatorUnit}} <: OperatorProd{V, I}

Operator.
"""
struct Operator{V<:Number, I<:ID{OperatorUnit}} <: OperatorProd{V, I}
    value::V
    id::I
end
@inline Operator(value::Number, id::OperatorUnit...) = Operator(value, id)
function Base.show(io::IO, m::Operator)
    @printf io "%s(%s%s%s)" nameof(typeof(m)) decimaltostr(value(m)) (rank(m)>0 ? ", " : "") join(id(m), ", ")
end


"""
    adjoint(m::Operator) -> Operator

Get the adjoint of an operator.
"""
@inline Base.adjoint(m::Operator) = rawtype(typeof(m))(value(m)', id(m)')



"""
    ishermitian(m::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
@inline ishermitian(m::Operator) = isa(value(m), Real) && ishermitian(id(m))


"""
    convert(::Type{M}, u::OperatorUnit) where {M<:Operator{<:Number, <:ID{OperatorUnit}}}

Convert an operator unit to an operator.
"""
@inline function Base.convert(::Type{M}, u::OperatorUnit) where {M<:Operator{<:Number, <:ID{OperatorUnit}}}
    @assert Tuple{typeof(u)} <: idtype(M) "convert error: not convertible."
    return Operator(one(valtype(M)), ID(u))
end








# Operator sum
"""
    OperatorSum{M<:OperatorPack, I<:Tuple} <: QuantumOperator

The sum of `OperatorPack`s.

Similar items are automatically merged with the aid of the id system.
"""
struct OperatorSum{M<:OperatorPack, I<:Tuple} <: QuantumOperator
    contents::Dict{I, M}
    OperatorSum(contents::Dict{<:Tuple, <:OperatorPack}) = new{valtype(contents), keytype(contents)}(contents)
end
@inline Base.eltype(ms::OperatorSum) = eltype(typeof(ms))
@inline Base.eltype(::Type{<:OperatorSum{M}}) where {M<:OperatorPack} = M
@inline Base.iterate(ms::OperatorSum) = iterate(values(ms.contents))
@inline Base.iterate(ms::OperatorSum, state) = iterate(values(ms.contents), state)
@inline Base.length(ms::OperatorSum) = length(ms.contents)
function Base.show(io::IO, ms::OperatorSum)
    @printf io "%s with %s %s\n" summary(ms) length(ms) nameof(eltype(ms))
    for m in ms
        @printf io "  %s\n" m
    end
end
@inline Base.haskey(ms::OperatorSum, id::Tuple) = haskey(ms.contents, id)
@inline Base.getindex(ms::OperatorSum, id::Tuple) = ms.contents[id]
@inline Base.setindex!(ms::OperatorSum, m::OperatorPack, id::Tuple) = (ms.contents[id] = m; m)
@inline Base.empty(ms::OperatorSum) = OperatorSum(empty(ms.contents))
@inline Base.empty!(ms::OperatorSum) = (empty!(ms.contents); ms)
@inline function Base.promote_rule(::Type{MS₁}, ::Type{MS₂}) where {MS₁<:OperatorSum, MS₂<:OperatorSum}
    M = promote_type(eltype(MS₁), eltype(MS₂))
    return OperatorSum{M, idtype(M)}
end
@inline function Base.convert(::Type{MS}, ms::OperatorSum) where {MS<:OperatorSum}
    @assert eltype(MS)>:eltype(ms) "convert error: cannot convert an object of $(typeof(ms)) to an object of $(MS)."
    return add!(zero(MS), ms)
end

"""
    OperatorSum(ms)
    OperatorSum(ms::QuantumOperator...)
    OperatorSum{M}(ms) where {M<:OperatorPack}
    OperatorSum{M}(ms::QuantumOperator...) where {M<:OperatorPack}

Get the sum of `OperatorPack`s.
"""
@inline OperatorSum(ms::QuantumOperator...) = OperatorSum{eltype(ms)}(ms)
@inline OperatorSum(ms) = OperatorSum{eltype(ms)}(ms)
@inline OperatorSum{M}(ms::QuantumOperator...) where {M<:OperatorPack} = OperatorSum{M}(ms)
function OperatorSum{M}(ms) where {M<:OperatorPack}
    result = OperatorSum(Dict{idtype(M), M}())
    for m in ms
        add!(result, m)
    end
    return result
end

"""
    isapprox(ms₁::OperatorSum, ms₂::OperatorSum; atol::Real=atol, rtol::Real=rtol) -> Bool

Compare two `OperatorSum`s and judge whether they are approximate to each other.
"""
@inline function Base.isapprox(ms₁::OperatorSum, ms₂::OperatorSum; atol::Real=atol, rtol::Real=rtol)
    for m in ms₁
        isapprox(value(m), 0, atol=atol, rtol=rtol) && continue
        k = id(m)
        haskey(ms₂, k) || return false
        isapprox(m, ms₂[k]) || return false
    end
    for m in ms₂
        isapprox(value(m), 0,  atol=atol, rtol=rtol) && continue
        k = id(m)
        haskey(ms₁, k) || return false
        isapprox(m, ms₁[k]) || return false
    end
    return true
end

"""
    zero(ms::OperatorSum) -> OperatorSum
    zero(::Type{MS}) where {MS<:OperatorSum} -> OperatorSum

Get the zero sum.
"""
@inline Base.zero(ms::OperatorSum) = zero(typeof(ms))
@inline Base.zero(::Type{MS}) where {MS<:OperatorSum} = OperatorSum{eltype(MS)}()

"""
    add!(ms::OperatorSum) -> typeof(ms)
    add!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack}) -> typeof(ms)
    add!(ms::OperatorSum, mms::OperatorSum) -> typeof(ms)

Get the in-place addition of quantum operators.
"""
@inline add!(ms::OperatorSum) = ms
@inline function add!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack})
    isa(m, Number) && m==0 && return ms
    m = convert(eltype(ms), m)
    old = get(ms.contents, id(m), nothing)
    new = isnothing(old) ? m : replace(old, value(old)+value(m))
    value(new)==0 ? delete!(ms.contents, id(m)) : (ms[id(m)] = new)
    return ms
end
@inline function add!(ms::OperatorSum, mms::OperatorSum)
    for m in mms
        add!(ms, m)
    end
    return ms
end

"""
    sub!(ms::OperatorSum) -> typeof(ms)
    sub!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack}) -> typeof(ms)
    sub!(ms::OperatorSum, mms::OperatorSum) -> typeof(ms)

Get the in-place subtraction of quantum operators.
"""
@inline sub!(ms::OperatorSum) = ms
@inline function sub!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack})
    isa(m, Number) && abs(m)==0 && return ms
    m = convert(eltype(ms), m)
    old = get(ms.contents, id(m), nothing)
    new = isnothing(old) ? -m : replace(old, value(old)-value(m))
    value(new)==0 ? delete!(ms.contents, id(m)) : (ms[id(m)] = new)
    return ms
end
@inline function sub!(ms::OperatorSum, mms::OperatorSum)
    for m in mms
        sub!(ms, m)
    end
    return ms
end

"""
    mul!(ms::OperatorSum, factor::Number) -> OperatorSum

Get the in-place multiplication of an `OperatorSum` with a number.
"""
function mul!(ms::OperatorSum, factor::Number)
    abs(factor)==0 && return empty!(ms)
    for m in ms
        new = replace(m, value(m)*factor)
        value(new)==0 ? delete!(ms.contents, id(m)) : (ms[id(m)] = new)
    end
    return ms
end

"""
    div!(ms::OperatorSum, factor::Number) -> OperatorSum

Get the in-place division of an `OperatorSum` with a number.
"""
@inline div!(ms::OperatorSum, factor::Number) = mul!(ms, one(dtype(eltype(ms)))/factor)

"""
    optype(m::QuantumOperator)
    optype(::Type{<:QuantumOperator})

Get the corresponding `OperatorPack` type of a generic quantum operator.
"""
@inline optype(m::QuantumOperator) = optype(typeof(m))
@inline optype(::Type{M}) where {M<:OperatorUnit} = fulltype(Operator, NamedTuple{(:value, :id), Tuple{Int, Tuple{M}}})
@inline optype(::Type{M}) where {M<:OperatorPack} = M
@inline optype(::Type{M}) where {M<:OperatorSum} = eltype(M)

"""
    +(m::QuantumOperator) -> typeof(m)
    +(m₁::QuantumOperator, m₂::QuantumOperator) -> OperatorSum
    +(factor::Number, m::QuantumOperator) -> OperatorSum
    +(m::QuantumOperator, factor::Number) -> OperatorSum

Overloaded `+` between quantum operators.
"""
@inline Base.:+(m::QuantumOperator) = m
@inline function Base.:+(m₁::QuantumOperator, m₂::QuantumOperator)
    M = promote_type(optype(m₁), optype(m₂))
    result = OperatorSum{M}()
    add!(result, m₁)
    add!(result, m₂)
    return result
end
@inline Base.:+(factor::Number, m::QuantumOperator) = one(optype(m))*factor + m
@inline Base.:+(m::QuantumOperator, factor::Number) = m + one(optype(m))*factor

"""
    -(m::QuantumOperator) -> QuantumOperator
    -(m₁::QuantumOperator, m₂::QuantumOperator) -> OperatorSum
    -(factor::Number, m::QuantumOperator) -> OperatorSum
    -(m::QuantumOperator, factor::Number) -> OperatorSum

Overloaded `-` between quantum operators.
"""
@inline Base.:-(m::QuantumOperator) = m*(-1)
@inline function Base.:-(m₁::QuantumOperator, m₂::QuantumOperator)
    M = promote_type(optype(m₁), optype(m₂))
    result = OperatorSum{M}()
    add!(result, m₁)
    sub!(result, m₂)
    return result
end
@inline Base.:-(factor::Number, m::QuantumOperator) = one(optype(m))*factor - m
@inline Base.:-(m::QuantumOperator, factor::Number) = m - one(optype(m))*factor

"""
    *(factor::Number, m::OperatorUnit) -> Operator
    *(m::OperatorUnit, factor::Number) -> Operator
    *(m₁::OperatorUnit, m₂::OperatorUnit) -> Operator
    *(factor::Number, m::OperatorPack) -> OperatorPack
    *(m::OperatorPack, factor::Number) -> OperatorPack
    *(m₁::OperatorPack, m₂::OperatorUnit) -> OperatorPack
    *(m₁::OperatorUnit, m₁::OperatorPack) -> OperatorPack
    *(m₁::OperatorPack, m₂::OperatorPack) -> OperatorPack
    *(factor::Number, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, factor::Number) -> OperatorSum
    *(m::OperatorPack, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, m::OperatorPack) -> OperatorSum
    *(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `*` between quantum operators or a quantum operator and a number.
"""
@inline Base.:*(factor::Number, m::OperatorUnit) = Operator(factor, m)
@inline Base.:*(m::OperatorUnit, factor::Number) = Operator(factor, m)
@inline Base.:*(m₁::OperatorUnit, m₂::OperatorUnit) = Operator(1, m₁, m₂)
@inline Base.:*(factor::Number, m::OperatorPack) = replace(m, factor*value(m))
@inline Base.:*(m::OperatorPack, factor::Number) = replace(m, value(m)*factor)
@inline Base.:*(m₁::OperatorPack, m₂::OperatorUnit) = m₁*Operator(1, m₂)
@inline Base.:*(m₁::OperatorUnit, m₂::OperatorPack) = Operator(1, m₁)*m₂
@inline function Base.:*(m₁::OperatorPack, m₂::OperatorPack)
    M₁, M₂ = typeof(m₁), typeof(m₂)
    @assert nameof(M₁)==nameof(M₂) && contentnames(M₁)==(:value, :id)==contentnames(M₂) "\"*\" error: not implemented between $(nameof(M₁)) and $(nameof(M₂))."
    return rawtype(M₁)(value(m₁)*value(m₂), ID(id(m₁), id(m₂)))
end
@inline Base.:*(factor::Number, ms::OperatorSum) = ms * factor
function Base.:*(ms::OperatorSum, factor::Number)
    abs(factor)==0 && return zero(ms)
    result = OperatorSum{promote_type(eltype(ms), typeof(factor))}()
    for m in ms
        add!(result, m*factor)
    end
    return result
end
@inline Base.:*(m::OperatorPack, ms::OperatorSum) = OperatorSum(collect(m*mm for mm in ms))
@inline Base.:*(ms::OperatorSum, m::OperatorPack) = OperatorSum(collect(mm*m for mm in ms))
@inline Base.:*(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum(collect(m₁*m₂ for m₁ in ms₁ for m₂ in ms₂))

"""
    /(m::QuantumOperator, factor::Number) -> QuantumOperator

Overloaded `/` between a quantum operator and a number.
"""
@inline Base.:/(m::QuantumOperator, factor::Number) = m * (one(dtype(optype(m)))/factor)

"""
    //(m::QuantumOperator, factor::Number) -> QuantumOperator

Overloaded `//` between a quantum operator and a number.
"""
@inline Base.://(m::QuantumOperator, factor::Number) = m * (1//factor)

"""
    ^(m::QuantumOperator, n::Integer) -> QuantumOperator

Overloaded `^` between a quantum operator and an integer.
"""
@inline Base.:^(m::QuantumOperator, n::Integer) = (@assert n>0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->m, Val(n)), init=1))

"""
    ⊗(m::Union{OperatorUnit, OperatorPack}, ms::OperatorSum) -> OperatorSum
    ⊗(ms::OperatorSum, m::Union{OperatorUnit, OperatorPack}) -> OperatorSum
    ⊗(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `⊗` between quantum operators.
"""
@inline ⊗(m::Union{OperatorUnit, OperatorPack}, ms::OperatorSum) = OperatorSum(collect(m ⊗ mm for mm in ms))
@inline ⊗(ms::OperatorSum, m::Union{OperatorUnit, OperatorPack}) = OperatorSum(collect(mm ⊗ m for mm in ms))
@inline ⊗(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum(collect(m₁ ⊗ m₂ for m₁ in ms₁ for m₂ in ms₂))

"""
    ⋅(m::Union{OperatorUnit, OperatorPack}, ms::OperatorSum) -> OperatorSum
    ⋅(ms::OperatorSum, m::Union{OperatorUnit, OperatorPack}) -> OperatorSum
    ⋅(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `⋅` between quantum operators.
"""
@inline ⋅(m::Union{OperatorUnit, OperatorPack}, ms::OperatorSum) = OperatorSum(collect(m ⋅ mm for mm in ms))
@inline ⋅(ms::OperatorSum, m::Union{OperatorUnit, OperatorPack}) = OperatorSum(collect(mm ⋅ m for mm in ms))
@inline ⋅(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum(collect(m₁ ⋅ m₂ for m₁ in ms₁ for m₂ in ms₂))










# Operators
"""
    Operators{O<:Operator, I<:ID{OperatorUnit}}

A set of operators.

Type alias for `OperatorSum{O<:Operator, I<:ID{OperatorUnit}}`.
"""
const Operators{O<:Operator, I<:ID{OperatorUnit}} = OperatorSum{O, I}
@inline Base.summary(io::IO, opts::Operators) = @printf io "Operators"

"""
    Operators(opts::Operator...)
    Operators{M}(opts::Operator...)

Get a set of operators.
"""
@inline Operators(opts::Operator...) = OperatorSum(opts)
@inline Operators{M}(opts::Operator...) where {M<:Operator} = OperatorSum{M}(opts)



"""
    adjoint(opts::Operators) -> Operators

Get the adjoint of a set of operators.
"""
function Base.adjoint(opts::Operators)
    result = zero(opts)
    for opt in opts
        add!(result, opt')
    end
    return result
end

"""
    ishermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
@inline ishermitian(opts::Operators) = opts == opts'








import LaTeXStrings: latexstring

"""
    latexstring(coupling::Coupling) -> String

Convert a `Coupling` to the latex format.
"""
function latexstring(coupling::Coupling)
    result = String[]
    len, count = length(coupling.constraint.representations), 1
    for i = 1:len
        r = rank(coupling.constraint, i)
        start, stop = count, count+r-1
        indexes = coupling.indexes[start:stop]
        pattern = coupling.constraint.representations[i]
        summation = pattern=="pattern" ? "" : replace(replace("$(coupling.constraint.representations[i])", "&&"=>"\\,\\text{and}\\,"), "||"=>"\\,\\text{or}\\,")
        summation=="" || (summation = join(push!(symbols(indexes, pattern), summation), ",\\,"))
        temp = isdefinite(indexes) ? "" : @sprintf "\\sum_{%s}" summation
        for index in indexes
            temp = @sprintf "%s %s" temp latexstring(index)
        end
        push!(result, temp)
        count = stop+1
    end
    return @sprintf "%s%s" valuetostr(coupling.value) join(result, " \\cdot ")
end

function symbols(indexes::Tuple{Vararg{Index}}, constraint::String)
    result = String[]
    for index in indexes
        iid = index.iid
        for i = 1:fieldcount(typeof(iid))
            attr = getfield(iid, i)
            if isa(attr, Symbol)
                value = string(attr)
                occursin(value, constraint) || push!(result, value)
            end
        end
    end
    return sort!(unique!(result))
end












# Term
"""
    TermFunction <: Function

Abstract type for concrete term functions.
"""
abstract type TermFunction <: Function end
@inline Base.:(==)(tf₁::TermFunction, tf₂::TermFunction) = ==(efficientoperations, tf₁, tf₂)
@inline Base.isequal(tf₁::TermFunction, tf₂::TermFunction) = isequal(efficientoperations, tf₁, tf₂)












"""
    Bond{K, P<:Point}

A generic bond, which could contains several points.
"""
struct Bond{K, P<:Point}
    kind::K
    points::Vector{P}
end
@inline Base.:(==)(bond₁::Bond, bond₂::Bond) = ==(efficientoperations, bond₁, bond₂)
@inline Base.isequal(bond₁::Bond, bond₂::Bond) = isequal(efficientoperations, bond₁, bond₂)
@inline Base.show(io::IO, bond::Bond) = @printf io "Bond(%s, %s)" repr(bond.kind) join(map(string, bond.points), ", ")

"""
    Bond(point::Point)
    Bond(kind, point₁::Point, point₂::Point, points::Point...)

Construct a bond.
"""
@inline Bond(point::Point) = Bond(0, [point])
@inline Bond(kind, point₁::Point, point₂::Point, points::Point...) = Bond(kind, [point₁, point₂, points...])










"""
    TermCoupling{E<:Coupling, C} <: TermFunction

The function for the coupling of a term.
"""
struct TermCoupling{E<:Coupling, C} <: TermFunction
    coupling::C
    TermCoupling(coupling) = new{eltype(coupling), typeof(coupling)}(coupling)
    TermCoupling(coupling::Function) = new{eltype(commontype(coupling, Tuple{Bond}, Any)), typeof(coupling)}(coupling)
    TermCoupling{E}(coupling::Function) where {E<:Coupling} = new{E, typeof(coupling)}(coupling)
end
@inline Base.valtype(termcoupling::TermCoupling) = valtype(typeof(termcoupling))
@inline Base.valtype(::Type{<:TermCoupling{E}}) where {E<:Coupling} = E
@inline (termcoupling::TermCoupling)(::Bond) = termcoupling.coupling
@inline (termcoupling::TermCoupling{<:Coupling, <:Function})(bond::Bond) = termcoupling.coupling(bond)






"""
    TermModulate(id::Symbol, modulate::Function)
    TermModulate(id::Symbol, modulate::Bool)

The function for the modulation of a term.
"""
struct TermModulate{M<:Union{Function, Val{true}, Val{false}}, id} <: TermFunction
    modulate::M
    TermModulate(id::Symbol, modulate::Function) = new{typeof(modulate), id}(modulate)
    TermModulate(id::Symbol, modulate::Bool=true) = new{Val{modulate}, id}(modulate|>Val)
end
@inline (termmodulate::TermModulate{Val{true}, id})(args...; kwargs...) where id = get(kwargs, id, nothing)
@inline (termmodulate::TermModulate{<:Function})(args...; kwargs...) = termmodulate.modulate(args...; kwargs...)
@inline ismodulatable(termmodulate::TermModulate) = ismodulatable(typeof(termmodulate))
@inline ismodulatable(::Type{<:TermModulate{Val{B}}}) where B = B
@inline ismodulatable(::Type{<:TermModulate{<:Function}}) = true










"""
    TermAmplitude(amplitude::Union{Function, Nothing}=nothing)

The function for the amplitude of a term.
"""
struct TermAmplitude{A<:Union{Function, Nothing}} <: TermFunction
    amplitude::A
    TermAmplitude(amplitude::Union{Function, Nothing}=nothing) = new{typeof(amplitude)}(amplitude)
end
@inline (termamplitude::TermAmplitude{Nothing})(::Bond) = 1
@inline (termamplitude::TermAmplitude{<:Function})(bond::Bond) = termamplitude.amplitude(bond)








"""
    Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}

A term of a quantum lattice system.
"""
mutable struct Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}
    value::V
    const bondkind::B
    const coupling::C
    const amplitude::A
    const ishermitian::Bool
    const modulate::M
    const factor::V
    function Term{K, I}(value, bondkind, coupling::TermCoupling, amplitude::TermAmplitude, ishermitian::Bool, modulate::TermModulate, factor) where {K, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        @assert value==value' "Term error: only real values are allowed. Complex values should be specified by the amplitude function."
        V, B, C, A, M = typeof(value), typeof(bondkind), typeof(coupling), typeof(amplitude), typeof(modulate)
        new{K, I, V, B, C, A, M}(value, bondkind, coupling, amplitude, ishermitian, modulate, factor)
    end
end
@inline Base.:(==)(term₁::Term, term₂::Term) = ==(efficientoperations, term₁, term₂)
@inline Base.isequal(term₁::Term, term₂::Term) = isequal(efficientoperations, term₁, term₂)

"""
    Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true) where K

Construct a term.
"""
@inline function Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true) where K
    return Term{K, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), ishermitian, TermModulate(id, modulate), 1)
end

"""
    kind(term::Term) -> Symbol
    kind(::Type{<:Term) -> Symbol

Get the kind of a term.
"""
@inline kind(term::Term) = kind(typeof(term))
@inline kind(::Type{<:Term{K}}) where K = K

"""
    id(term::Term) -> Symbol
    id(::Type{<:Term) -> Symbol

Get the id of a term.
"""
@inline id(term::Term) = id(typeof(term))
@inline id(::Type{<:Term{K, I} where K}) where I = I

"""
    valtype(term::Term)
    valtype(::Type{<:Term)

Get the value type of a term.
"""
@inline Base.valtype(term::Term) = valtype(typeof(term))
@inline Base.valtype(::Type{<:Term{K, I, V} where {K, I}}) where V = V

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term) -> Int

Get the rank of a term.
"""
@inline rank(term::Term) = rank(typeof(term))
@inline rank(::Type{<:Term{K, I, V, B, C} where {K, I, V, B}}) where {C<:TermCoupling} = rank(valtype(C))

"""
    ismodulatable(term::Term) -> Bool
    ismodulatable(::Type{<:Term}) -> Bool

Judge whether a term could be modulated by its modulate function.
"""
@inline ismodulatable(term::Term) = ismodulatable(typeof(term))
@inline ismodulatable(::Type{<:Term{K, I, V, B, <:TermCoupling, <:TermAmplitude, M} where {K, I, V, B}}) where M = ismodulatable(M)

"""
    repr(term::Term, bond::Bond, hilbert::Hilbert) -> String

Get the repr representation of a term on a bond with a given Hilbert space.
"""
function Base.repr(term::Term, bond::Bond, hilbert::Hilbert)
    cache = String[]
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                if !isnothing(iterate(expand(term|>kind|>Val, coupling, bond, hilbert)))
                    representation = repr(value*coupling)
                    term.ishermitian || (representation = string(replace(representation, " * "=>"*"), " + h.c."))
                    push!(cache, @sprintf "%s" representation)
                end
            end
        end
    end
    return join(cache, "\n")
end

"""
    replace(term::Term; kwargs...) -> Term

Replace some attributes of a term with key word arguments.
"""
@inline @generated function Base.replace(term::Term; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(term, $name))) for name in QuoteNode.(term|>fieldnames)]
    return :(Term{kind(term), id(term)}($(exprs...)))
end











"""
    Onsite(id::Symbol, value, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Onsite term.

Type alias for `Term{:Onsite, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Onsite{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Onsite, id, V, Int, C, A, M}
@inline function Onsite(id::Symbol, value, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Onsite}(id, value, 0, coupling, ishermitian; amplitude=amplitude, modulate=modulate)
end














"""
    Hopping(id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Hopping term.

Type alias for `Term{:Hopping, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hopping{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hopping, id, V, B, C, A, M}
@inline function Hopping(id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    return Term{:Hopping}(id, value, bondkind, coupling, false; amplitude=amplitude, modulate=modulate)
end





"""
    Pairing(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Pairing term.

Type alias for `Term{:Pairing, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Pairing{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Pairing, id, V, B, C, A, M}
@inline function Pairing(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Pairing}(id, value, bondkind, coupling, false; amplitude=amplitude, modulate=modulate)
end
@inline nambu(::Val{:Pairing}, ::Colon, ::Int) = annihilation
function expand!(operators::Operators, term::Pairing, bond::Bond, hilbert::Hilbert; half::Bool=false)
    argtypes = Tuple{Operators, Term, Bond, Hilbert}
    invoke(expand!, argtypes, operators, term, bond, hilbert; half=half)
    length(bond)==2 && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half)
    return operators
end








"""
    Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Hubbard term.

Type alias for `Term{:Hubbard, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hubbard{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hubbard, id, V, Int, C, A, M}
@inline function Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Hubbard}(id, value, 0, Coupling(:, FID, :, (1//2, 1//2, -1//2, -1//2), (2, 1, 2, 1)), true; amplitude=amplitude, modulate=modulate)
end






"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalInterSpin, id, V, Int, C, A, M}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:InterOrbitalInterSpin}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, σ₁, 2)), Index(:, FID(α, σ₁, 1)), Index(:, FID(β, σ₂, 2)), Index(:, FID(β, σ₂, 1)); constraint=α<β && σ₁≠σ₂)), true;
        amplitude=amplitude, modulate=modulate
    )
end





"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalIntraSpin, id, V, Int, C, A, M}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:InterOrbitalIntraSpin}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, σ, 2)), Index(:, FID(α, σ, 1)), Index(:, FID(β, σ, 2)), Index(:, FID(β, σ, 1)); constraint=α<β)), true;
        amplitude=amplitude, modulate=modulate
    )
end










"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Spin-flip term.

Type alias for `Term{:SpinFlip, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinFlip{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:SpinFlip, id, V, Int, C, A, M}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:SpinFlip}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, 1//2, 2)), Index(:, FID(β, -1//2, 2)), Index(:, FID(α, -1//2, 1)), Index(:, FID(β, 1//2, 1)); constraint=α<β)), false;
        amplitude=amplitude, modulate=modulate
    )
end








"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Pair-hopping term.

Type alias for `Term{:PairHopping, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const PairHopping{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:PairHopping, id, V, Int, C, A, M}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:PairHopping}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, 1//2, 2)), Index(:, FID(α, -1//2, 2)), Index(:, FID(β, -1//2, 1)), Index(:, FID(β, 1//2, 1)); constraint=α<β)), false;
        amplitude=amplitude, modulate=modulate
    )
end








"""
    Coulomb(
        id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :)))^2;
        ishermitian::Bool=true,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=true
    )

Coulomb term.

Type alias for `Term{:Coulomb, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Coulomb{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Coulomb, id, V, B, C, A, M}
@inline function Coulomb(
    id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :)))^2;
    ishermitian::Bool=true,
    amplitude::Union{Function, Nothing}=nothing,
    modulate::Union{Function, Bool}=true
)
    return Term{:Coulomb}(id, value, bondkind, coupling, ishermitian; amplitude=amplitude, modulate=modulate)
end







"""
    FockTerm

Type alias for `Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}`.
"""
const FockTerm = Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}









"""
    Parameters{Names}(values::Number...) where Names

A NamedTuple that contain the key-value pairs.
"""
const Parameters{Names} = NamedTuple{Names, <:Tuple{Vararg{Number}}}
@inline Parameters{Names}(values::Number...) where {Names} = NamedTuple{Names}(values)
@inline @generated function update(params::NamedTuple; parameters...)
    names = fieldnames(params)
    values = Expr(:tuple, [:(get(parameters, $name, getfield(params, $name))) for name in QuoteNode.(names)]...)
    return :(NamedTuple{$names}($values))
end














# Transformations
"""
    Transformation <: Function

Abstract transformation on quantum operators.
"""
abstract type Transformation <: Function end
@inline Base.:(==)(transformation₁::Transformation, transformation₂::Transformation) = ==(efficientoperations, transformation₁, transformation₂)
@inline Base.isequal(transformation₁::Transformation, transformation₂::Transformation) = isequal(efficientoperations, transformation₁, transformation₂)

"""
    LinearTransformation <: Transformation

Abstract linear transformation on quantum operators.
"""
abstract type LinearTransformation <: Transformation end
@inline Base.valtype(transformation::LinearTransformation, m::QuantumOperator) = valtype(typeof(transformation), typeof(m))
@inline Base.zero(transformation::LinearTransformation, m::QuantumOperator) = zero(valtype(transformation, m))











# Boundary
"""
    Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names

Boundary twist of operators.
"""
struct Boundary{Names, D<:Number, V<:AbstractVector} <: LinearTransformation
    values::Vector{D}
    vectors::Vector{V}
    function Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names
        @assert length(Names)==length(values)==length(vectors) "Boundary error: mismatched names, values and vectors."
        datatype = promote_type(eltype(values), Float)
        new{Names, datatype, eltype(vectors)}(convert(Vector{datatype}, values), vectors)
    end
end
@inline Base.:(==)(bound₁::Boundary, bound₂::Boundary) = keys(bound₁)==keys(bound₂) && ==(efficientoperations, bound₁, bound₂)
@inline Base.isequal(bound₁::Boundary, bound₂::Boundary) = isequal(keys(bound₁), keys(bound₂)) && isequal(efficientoperations, bound₁, bound₂)
@inline Base.valtype(::Type{<:Boundary}, M::Type{<:Operator}) = reparameter(M, :value, promote_type(Complex{Int}, valtype(M)))
@inline Base.valtype(B::Type{<:Boundary}, MS::Type{<:Operators}) = (M = valtype(B, eltype(MS)); Operators{M, idtype(M)})

"""
    keys(bound::Boundary) -> Tuple{Vararg{Symbol}}
    keys(::Type{<:Boundary{Names}}) where Names -> Names

Get the names of the boundary parameters.
"""
@inline Base.keys(bound::Boundary) = keys(typeof(bound))
@inline Base.keys(::Type{<:Boundary{Names}}) where Names = Names






"""
    Frontend

The frontend of algorithms applied to a quantum lattice system.
"""
abstract type Frontend end
@inline Base.:(==)(frontend₁::Frontend, frontend₂::Frontend) = ==(efficientoperations, frontend₁, frontend₂)
@inline Base.isequal(frontend₁::Frontend, frontend₂::Frontend) = isequal(efficientoperations, frontend₁, frontend₂)
@inline Base.repr(frontend::Frontend) = String(nameof(typeof(frontend)))
@inline Base.show(io::IO, frontend::Frontend) = @printf io "%s" nameof(typeof(frontend))
@inline Base.valtype(frontend::Frontend) = valtype(typeof(frontend))
@inline update!(frontend::Frontend; kwargs...) = error("update! error: not implemented for $(nameof(typeof(frontend))).")
@inline Parameters(frontend::Frontend) = error("Parameters error: not implemented for $(nameof(typeof(frontend))).")







"""
    RepresentationGenerator <: Frontend

Representation generator of a quantum lattice system.
"""
abstract type RepresentationGenerator <: Frontend end
@inline Base.eltype(gen::RepresentationGenerator) = eltype(typeof(gen))
@inline Base.IteratorSize(::Type{<:RepresentationGenerator}) = Base.SizeUnknown()
@inline function Base.iterate(gen::RepresentationGenerator)
    ops = expand(gen)
    index = iterate(ops)
    isnothing(index) && return nothing
    return index[1], (ops, index[2])
end
@inline function Base.iterate(gen::RepresentationGenerator, state)
    index = iterate(state[1], state[2])
    isnothing(index) && return nothing
    return index[1], (state[1], index[2])
end
@inline Base.show(io::IO, ::MIME"text/latex", gen::RepresentationGenerator) = show(io, MIME"text/latex"(), latexstring(latexstring(expand(gen))))

"""
    expand(gen::RepresentationGenerator) -> valtype(gen)
    expand!(result, gen::RepresentationGenerator) -> typeof(result)

Expand the generator to get the representation of the quantum lattice system (or some part of it).
"""
@inline expand(gen::RepresentationGenerator) = expand!(zero(valtype(gen)), gen)
@inline expand!(result, gen::RepresentationGenerator) = error("expand! error: not implemented for $(nameof(typeof(gen))).")
















"""
    Entry{C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary} <: RepresentationGenerator

The basic representation generator of a quantum lattice system that records the quantum operators or a representation of the quantum operators related to (part of) the system.
"""
mutable struct Entry{C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary} <: RepresentationGenerator
    const constops::C
    const alterops::A
    const boundops::B
    parameters::P
    const boundary::D
end
@inline Entry(entry::Entry) = entry
@inline Base.eltype(E::Type{<:Entry}) = eltype(valtype(E))




"""
    Entry(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain)

Construct an entry of quantum operators based on the input terms, bonds, Hilbert space and (twisted) boundary condition.
"""
function Entry(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain)
    constterms, alterterms, choosedterms = termclassifier(terms)
    innerbonds = filter(bond->isintracell(bond), bonds)
    boundbonds = filter(bond->!isintracell(bond), bonds)
    constops = Operators{mapreduce(term->optype(typeof(term), typeof(hilbert), eltype(bonds)), promote_type, choosedterms)}()
    map(term->expand!(constops, term, innerbonds, hilbert; half=half), constterms)
    alterops = NamedTuple{map(id, alterterms)}(map(term->expand(one(term), innerbonds, hilbert; half=half), alterterms))
    boundops = NamedTuple{map(id, terms)}(map(term->map!(boundary, expand!(Operators{valtype(typeof(boundary), optype(typeof(term), typeof(hilbert), eltype(bonds)))}(), one(term), boundbonds, hilbert, half=half)), terms))
    parameters = NamedTuple{map(id, terms)}(map(term->term.value, terms))
    return Entry(constops, alterops, boundops, parameters, boundary)
end

@generated function termclassifier(terms::Tuple{Vararg{Term}})
    constterms, alterterms = [], []
    for (i, term) in enumerate(fieldtypes(terms))
        ismodulatable(term) ? push!(alterterms, :(terms[$i])) : push!(constterms, :(terms[$i]))
    end
    constterms, alterterms = Expr(:tuple, constterms...), Expr(:tuple, alterterms...)
    return Expr(:tuple, constterms, alterterms, (length(constterms.args)>0 ? constterms : alterterms))
end







# Metric and Table
"""
    Metric <: Function

The rules for measuring an operator unit so that different operator units can be compared.

As a function, every instance should accept only one positional argument, i.e. the operator unit to be measured.
"""
abstract type Metric <: Function end
@inline Base.:(==)(m₁::T, m₂::T) where {T<:Metric} = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::T, m₂::T) where {T<:Metric} = isequal(efficientoperations, m₁, m₂)
@inline (M::Type{<:Metric})(::Type{I}) where {I<:AbstractCompositeIndex} = M(indextype(I))
@inline (metric::Metric)(index::AbstractCompositeIndex) = metric(getcontent(index, :index))
@inline Base.valtype(::Type{M}, ::Type{I}) where {M<:Metric, I<:AbstractCompositeIndex} = valtype(M, indextype(I))
@inline (M::Type{<:Metric})(::Type{H}) where {H<:Hilbert} = M(Index{Int, H|>valtype|>eltype})

"""
    OperatorUnitToTuple{Fields} <: Metric

A rule that converts an operator unit to a tuple by iterating over a set of selected fields in a specific order.
"""
struct OperatorUnitToTuple{Fields} <: Metric
    OperatorUnitToTuple(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline OperatorUnitToTuple(fields::Symbol...) = OperatorUnitToTuple(fields)

"""
    keys(::OperatorUnitToTuple{Fields}) where Fields -> Fields
    keys(::Type{<:OperatorUnitToTuple{Fields}}) where Fields -> Fields

Get the names of the selected fields.
"""
@inline Base.keys(::OperatorUnitToTuple{Fields}) where Fields = Fields
@inline Base.keys(::Type{<:OperatorUnitToTuple{Fields}}) where Fields = Fields

"""
    OperatorUnitToTuple(::Type{I}) where {I<:Index}

Construct the metric rule from the information of the `Index` type.
"""
@inline OperatorUnitToTuple(::Type{I}) where {I<:Index} = OperatorUnitToTuple(:site, (fieldnames(iidtype(I)))...)

"""
    valtype(::Type{<:OperatorUnitToTuple}, ::Type{<:Index})

Get the valtype of applying an `OperatorUnitToTuple` rule to an `Index`.
"""
@inline @generated function Base.valtype(::Type{M}, ::Type{I}) where {M<:OperatorUnitToTuple, I<:Index}
    types = []
    for field in keys(M)
        if field==:site
            push!(types, Int)
        elseif hasfield(iidtype(I), field)
            push!(types, fieldtype(iidtype(I), field))
        end
    end
    return  Expr(:curly, :Tuple, types...)
end

"""
    (operatorunittotuple::OperatorUnitToTuple)(index::Index) -> Tuple

Convert an index to a tuple.
"""
@inline @generated function (operatorunittotuple::OperatorUnitToTuple)(index::Index)
    exprs = []
    for name in keys(operatorunittotuple)
        field = QuoteNode(name)
        if name==:site
            push!(exprs, :(index.site))
        elseif hasfield(iidtype(index), name)
            push!(exprs, :(getfield(index.iid, $field)))
        end
    end
    return Expr(:tuple, exprs...)
end




"""
    Table{I, B<:Metric} <: CompositeDict{I, Int}

The table of operator unit vs. sequence pairs.
"""
struct Table{I, B<:Metric} <: CompositeDict{I, Int}
    by::B
    contents::Dict{I, Int}
end
@inline contentnames(::Type{<:Table}) = (:by, :contents)
@inline Table{I}(by::Metric) where {I<:OperatorUnit} = Table(by, Dict{valtype(typeof(by), I), Int}())
@inline vec2dict(vs::AbstractVector) = Dict{eltype(vs), Int}(v=>i for (i, v) in enumerate(vs))

"""
    getindex(table::Table, operatorunit::OperatorUnit) -> Int

Inquiry the sequence of an operator unit.
"""
@inline Base.getindex(table::Table, operatorunit::OperatorUnit) = table[table.by(operatorunit)]

"""
    haskey(table::Table, operatorunit::OperatorUnit) -> Bool
    haskey(table::Table, operatorunits::ID{OperatorUnit}) -> Tuple{Vararg{Bool}}

Judge whether a single operator unit or a set of operator units have been assigned with sequences in table.
"""
@inline Base.haskey(table::Table, operatorunit::OperatorUnit) = haskey(table, table.by(operatorunit))
@inline Base.haskey(table::Table, operatorunits::ID{OperatorUnit}) = map(operatorunit->haskey(table, operatorunit), operatorunits)

"""
    Table(operatorunits::AbstractVector{<:OperatorUnit}, by::Metric=OperatorUnitToTuple(eltype(operatorunits)))

Convert a set of operator units to the corresponding table of operator unit vs. sequence pairs.

The input operator units are measured by the input `by` function with the duplicates removed. The resulting unique 
values are sorted, which determines the sequence of the input `operatorunits`. Note that two operator units have 
the same sequence if their converted values are equal to each other.
"""
@inline Table(operatorunits::AbstractVector{<:OperatorUnit}, by::Metric=OperatorUnitToTuple(eltype(operatorunits))) = Table(by, [by(operatorunit) for operatorunit in operatorunits]|>unique!|>sort!|>vec2dict)

"""
    Table(hilbert::Hilbert, by::Metric=OperatorUnitToTuple(typeof(hilbert))) -> Table

Get the index-sequence table of a Hilbert space.
"""
function Table(hilbert::Hilbert, by::Metric=OperatorUnitToTuple(typeof(hilbert)))
    result = Index{Int, hilbert|>valtype|>eltype}[]
    for (site, internal) in hilbert
        for iid in internal
            push!(result, (result|>eltype)(site, iid))
        end
    end
    return Table(result, by)
end









"""
    OperatorGenerator{E<:Entry{<:Operators}, TS<:Tuple{Vararg{Term}}, B<:Bond, H<:Hilbert, T<:Union{Table, Nothing}} <: CompositeGenerator{E, T}

A generator of operators based on the terms, bonds and Hilbert space of a quantum lattice system.
"""
struct OperatorGenerator{E<:Entry{<:Operators}, TS<:Tuple{Vararg{Term}}, B<:Bond, H<:Hilbert, T<:Union{Table, Nothing}} <: CompositeGenerator{E, T}
    operators::E
    terms::TS
    bonds::Vector{B}
    hilbert::H
    half::Bool
    table::T
end
@inline contentnames(::Type{<:OperatorGenerator}) = (:operators, :terms, :bonds, :hilbert, :half, :table)



"""
    OperatorGenerator(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain, table::Union{Table,Nothing}=nothing)

Construct a generator of operators.
"""
@inline function OperatorGenerator(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain, table::Union{Table,Nothing}=nothing)
    return OperatorGenerator(Entry(terms, bonds, hilbert; half=half, boundary=boundary), terms, bonds, hilbert, half, table)
end





