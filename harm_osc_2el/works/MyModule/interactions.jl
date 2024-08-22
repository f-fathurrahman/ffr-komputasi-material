abstract type Hamiltonian end # need this??
abstract type Interaction <: Hamiltonian end


# Yukawa potential in 1d?
struct HarmonicOscillatorCoulomb <: Interaction
    ω2::Float64
    shielding2::Float64
    # Inner constructor, arguments need to be squared
    function HarmonicOscillatorCoulomb(ω, shielding)
        new(ω^2, shielding^2)
    end
end

# Need this???
HarmonicOscillatorCoulomb(ω; shielding=0.0) = HarmonicOscillatorCoulomb(ω, shielding)
# XXX: Default params for shielding

function evaluate_on_grid!(Vx1, x1, grid, V::HarmonicOscillatorCoulomb)
    Vx1 .= 1 ./ sqrt.( (grid .- x1).^2 .+ V.shielding2 )
    return
end