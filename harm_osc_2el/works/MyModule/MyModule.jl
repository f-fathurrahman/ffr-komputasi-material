module MyModule

using Infiltrator
using LinearAlgebra: Diagonal, kron

include("HarmonicOscillatorBasis.jl");
export HarmonicOscillatorBasis
export SpinBasis
export evaluate_basis
export evaluate_on_grid

include("interactions.jl")
export HarmonicOscillatorCoulomb

include("functions.jl")
export calc_onebody_integrals
export calc_twobody_integrals

end