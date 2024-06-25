using LinearAlgebra
using Formatting
import Base: length, eltype, show

export MyFermiSparse

abstract type AbstractMyFermiArray{T,N} <: AbstractArray{T,N} end

struct MyFermiSparse{Td,Ti,N} <: AbstractMyFermiArray{Td,N}
    indexes::Vector{NTuple{N, Ti}}
    data::Vector{Td}
end

function ndims(A::MyFermiSparse{Td,Ti,N}) where {Td, Ti, N}
    return N
end

show(io::IO, ::MIME"text/plain", A::MyFermiSparse{Td, Ti, N}) where {Td, Ti, N} = show(io,A)

function show(io::IO, A::MyFermiSparse{Td, Ti, N}) where {Td, Ti, N}
    println(io, "MyFermiSparse{$Td, $Ti, $N}:\n")
    str_out = ""
    for (idx, val) in zip(A.indexes, A.data)
        str_idx = "("
        for k in idx
            str_idx *= "$k,"
        end
        str_val = format("{:< 12.10f}", val)
        str_out *= str_idx[1:end-1]*") => $str_val\n"
    end
    print(io, str_out)
end

function index2(i::Signed, j::Signed)::Signed
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

# Basic methods for AbstractArrays in Julia
function length(A::AbstractMyFermiArray)
    return length(A.data)
end

function eltype(A::AbstractMyFermiArray)
    return eltype(A.data)
end