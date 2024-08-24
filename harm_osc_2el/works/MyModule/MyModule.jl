module MyModule

using Infiltrator
using Printf
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

include("rhf.jl")
export RHFState, RHF

include("calc_ground_state.jl")
export calc_ground_state!

end # module

