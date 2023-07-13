module ACEbase

using Reexport, StaticArrays

include("def.jl")

include("objectpools.jl")
import ACEbase.ObjectPools: acquire!, release!

abstract type AbstractBasis end
abstract type ACEBasis <: AbstractBasis end
abstract type OneParticleBasis{T} <: ACEBasis end
abstract type Discrete1pBasis{T} <: OneParticleBasis{T} end 
abstract type ScalarACEBasis <: ACEBasis end

abstract type AbstractState end

abstract type AbstractConfiguration end

abstract type AbstractContinuousState <: AbstractState end

abstract type AbstractDiscreteState <: AbstractState end

isdiscrete(::AbstractContinuousState) = false
isdiscrete(::AbstractDiscreteState) = true


"""
`function fltype`

Return the output floating point type employed by some object, typically a
calculator or basis.

DEPRECATE THIS!!!
"""
function fltype end

fltype(T::DataType) = T

fltype_intersect(o1, o2) =
   fltype_intersect(fltype(o1), fltype(o2))

fltype_intersect(T1::DataType, T2::DataType) =
   typeof(one(T1) * one(T2))

"""
`function rfltype end`

Return the real floating point type employed by some object,
typically a calculator or basis, this is normally the same as fltype, but
it can be difference e.g. `rfltype = real ∘ flype

DEPRECATE THIS!! 
"""
rfltype(args...) = real(fltype(args...))


# TODO: deprecate all of the following: 


# function alloc_B end

# alloc_B(basis::ACEBasis, x) = zeros(valtype(basis, x), length(basis))

# function alloc_dB end

# alloc_dB(B::ACEBasis, args...) = zeros(gradtype(B, args...), length(B))

# alloc_dB(B::ACEBasis, cfg::AbstractConfiguration) = 
#             zeros( gradtype(B, zero(eltype(cfg))), 
#                    (length(B), length(cfg) ) )

# -----------




function valtype end 

valtype(basis::ACEBasis, cfg::AbstractConfiguration) =
      valtype(basis, zero(eltype(cfg)))


function gradtype end




function combine end

function evaluate end
function evaluate_d end
function evaluate_dd end
function evaluate_ed end
function evaluate! end
function evaluate_d! end
function evaluate_dd! end
function evaluate_ed! end
function precon! end


evaluate(basis::ACEBasis, args...) =  
      evaluate!( acquire_B!(basis, args...), basis, args... )

evaluate_d(basis::ACEBasis, args...) =  
      evaluate_d!( acquire_dB!(basis, args...), basis, args... )

evaluate_ed(basis::ACEBasis, args...) =  
      evaluate_ed!( acquire_B!(basis, args...), acquire_dB!(basis, args...), 
                    basis, args... )



# TODO - documentation 
function acquire_B! end 
function release_B! end 
function acquire_dB! end 
function release_dB! end 

# a few nice fallbacks / defaults 

# function acquire_B!(basis::ACEBasis, args...) 
#    VT = valtype(basis, args...)
#    if hasproperty(basis, :B_pool)
#       return acquire!(basis.B_pool, length(basis), VT)
#    end
#    return Vector{VT}(undef, length(basis))
# end

# function release_B!(basis::ACEBasis, B) 
#    if hasproperty(basis, :B_pool)
#       release!(basis.B_pool, B)
#    end
# end 

@generated function acquire_B!(basis::TB, args...)  where {TB <: ACEBasis}
   if :B_pool in fieldnames(TB)
      return quote 
         VT = valtype(basis, args...)
         return acquire!(basis.B_pool, length(basis), VT)
      end
   else
      return quote 
         VT = valtype(basis, args...)
         Vector{VT}(undef, length(basis))::Vector{VT}
      end
   end
end

@generated function release_B!(basis::TB, B)  where {TB <: ACEBasis}
   if :B_pool in fieldnames(basis) 
      return quote 
         release!(basis.B_pool, B)
      end
   else
      return quote nothing; end 
   end
end 


# function acquire_dB!(basis::ACEBasis, args...) 
#    GT = gradtype(basis, args...)
#    if hasproperty(basis, :dB_pool)
#       return acquire!(basis.dB_pool, length(basis), GT)
#    end
#    return Vector{GT}(undef, length(basis))
# end

# function acquire_dB!(basis::ACEBasis, cfg::AbstractConfiguration) 
#    GT = gradtype(basis, cfg)
#    sz = (length(basis), length(cfg))
#    if hasproperty(basis, :dB_pool)
#       return acquire!(basis.dB_pool, sz, GT)
#    end
#    return Matrix{GT}(undef, sz)
# end

# function release_dB!(basis::ACEBasis, dB) 
#    if hasproperty(basis, :dB_pool)
#       release!(basis.dB_pool, dB)
#    end
# end 

@generated function acquire_dB!(basis::TB, args...)  where {TB <: ACEBasis}
   if :dB_pool in fieldnames(TB)
      return quote 
         GT = gradtype(basis, args...)
         acquire!(basis.dB_pool, length(basis), GT)
      end
   end
   return quote 
      GT = gradtype(basis, args...)
      return Vector{GT}(undef, length(basis))   
   end 
end

@generated function acquire_dB!(basis::TB, cfg::AbstractConfiguration)  where {TB <: ACEBasis}
   if :dB_pool in fieldnames(TB)
      return quote 
         GT = gradtype(basis, cfg)
         sz = (length(basis), length(cfg))
         return acquire!(basis.dB_pool, sz, GT)
      end
   end
   return quote 
      GT = gradtype(basis, cfg)
      sz = (length(basis), length(cfg))
      return Matrix{GT}(undef, sz)
   end 
end

@generated function release_dB!(basis::TB, B)  where {TB <: ACEBasis}
   if :dB_pool in fieldnames(basis)
      return quote 
         release!(basis.dB_pool, B)
      end
   else
      return quote nothing; end 
   end
end



# TODO: create a version where we can list which fields should be equal
"""
a simple utility function to check whether two objects are equal
"""
_allfieldsequal(x1, x2) =
      all( getfield(x1, sym) == getfield(x2, sym)
           for sym in union(fieldnames(typeof(x1)), fieldnames(typeof(x2))) )


include("fio.jl")
@reexport using ACEbase.FIO

include("testing.jl")


# ---- some utility code for discrete basis sets to avoid any attempt to 
#      differentiate it

# TODO: blech - nasty hack...
alloc_dB(basis::Discrete1pBasis{T}, N::Integer) where {T} = 
            zeros(SVector{3, T}, length(basis), N)





end
