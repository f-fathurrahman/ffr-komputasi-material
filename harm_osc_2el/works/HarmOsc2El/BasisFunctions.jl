abstract type Basis end

abstract type SpatialBasis <: Basis end

include("HObasis.jl")
include("spinbasis.jl")
include("pairingbasis.jl")
