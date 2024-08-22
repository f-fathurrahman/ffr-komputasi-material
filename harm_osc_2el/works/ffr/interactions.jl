abstract type Hamiltonian end # need this??
abstract type Interaction <: Hamiltonian end

struct HarmonicOscillatorCoulomb <: Interaction
    ω2::Float64
    shielding2::Float64
end
function HarmonicOscillatorCoulomb(ω, shielding)
    HarmonicOscillatorCoulomb(ω^2, shielding^2)
end
# Need this???
HarmonicOscillatorCoulomb(ω; shielding=0.0) = HarmonicOscillatorCoulomb(ω, shielding)
# XXX: Default params for shielding

function evaluate_on_grid!(interaction, x1, grid, V::HarmonicOscillatorCoulomb)
    interaction .= 1 ./ sqrt.( (grid .- x1).^2 .+ V.shielding2 )
    return interaction
end