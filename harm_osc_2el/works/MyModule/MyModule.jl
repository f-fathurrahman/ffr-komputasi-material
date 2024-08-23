module MyModule

using Infiltrator
import LinearAlgebra
using LinearAlgebra: Diagonal, kron, kron!, dot, pinv, diag

include("HarmonicOscillatorBasis.jl");
export SpatialBasis
export HarmonicOscillatorBasis
export SpinBasis
export evaluate_basis
export evaluate_on_grid

include("interactions.jl")
export HarmonicOscillatorCoulomb

include("integrals.jl")
export calc_onebody_integrals
export calc_twobody_integrals

include("system.jl")
export SpatialSystem, init_system

include("mixer.jl")
export LinearMixer, DIISMixer, compute_new_vector

end # module

