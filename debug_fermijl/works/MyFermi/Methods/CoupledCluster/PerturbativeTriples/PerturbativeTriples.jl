abstract type RpTAlgorithm end

function get_rpt_alg()
    implemented = [ijk(), ijk2(), abc()]
    N = Options.get("pt_alg")
    try 
        return implemented[N]
    catch BoundsError
        throw(MyFermiException("implementation number $N not available for RCCSD(T)."))
    end
end
"""
    MyFermi.CoupledCluster.RCCSDpT

MyFermi struct that holds information about RCCSD(T) wavefunctions

_struct tree:_

**RCCSD(T)** <: AbstractCCWavefunction <: AbstractCorrelatedWavefunction <: AbstractWavefunction
"""
struct RCCSDpT{T} <: AbstractCCWavefunction
    CCSD::RCCSD{T}
    energy::T
    correction::T
end

"""
    MyFermi.CoupledCluster.get_pt_alg

Returns a singleton type corresponding to a RCCSD(T) implementation based on the options.
"""
function get_pt_alg(N::Int = Options.get("pt_alg"))
    try 
        return get_pt_alg(Val(N))
    catch MethodError
        throw(MyFermiException("implementation number $N not available for RCCSD(T)."))
    end
end

# Implementations
struct ijk <: RpTAlgorithm end
include("ijk.jl")
get_pt_alg(x::Val{1}) = ijk()
struct ijk2 <: RpTAlgorithm end
include("ijk2.jl")
get_pt_alg(x::Val{2}) = ijk2()
struct abc <: RpTAlgorithm end
include("abc.jl")
get_pt_alg(x::Val{3}) = abc()

function RCCSDpT(x...)
    if !any(i-> i isa RpTAlgorithm, x)
        RCCSDpT(x..., get_rpt_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(MyFermiException("invalid arguments for RCCSD(T) method: $args"))
    end
end

# Gradient methods
include("Gradients/RCCSDpTgrad.jl")