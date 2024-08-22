module HarmOsc2El

import LinearAlgebra
using TensorOperations: @tensor

# This will include various other files
include("system.jl")

# Orbital systems
export System, Basis, reference_energy, particle_density
export SpatialSystem, SpatialBasis, HOBasis, SpinBasis
export PairingSystem, pairing_exact, pairing_MBPT2, Pairing
export CalogeroSutherland, HOCoulomb, HarmonicOscillator
export RHF, HF, CCD, CCSD
export energy
export compute_ground_state!

include("mixer.jl")
include("util_methods.jl")
include("hf.jl")
include("rhf.jl")
include("ccd.jl")
include("ccsd.jl")

end