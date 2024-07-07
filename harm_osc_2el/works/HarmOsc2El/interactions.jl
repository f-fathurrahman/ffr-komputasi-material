abstract type Hamiltonian end

abstract type Interaction <: Hamiltonian end

struct HOCoulomb <: Interaction
    ω2::Float64
    shielding2::Float64
    function HOCoulomb(ω, shielding)
        new(ω^2, shielding^2)
    end
end
HOCoulomb(ω; shielding) = HOCoulomb(ω, shielding)

function interaction_over_grid!(interaction, x1, grid, V::HOCoulomb)
    interaction .= 1 ./ sqrt.( (grid .- x1).^2 .+ V.shielding2 )
    return interaction
end

struct CalogeroSutherland <: Interaction
    ββ_1::Float64
    function CalogeroSutherland(beta) # 2 is normal
        ββ_1 = beta * (beta - 1)
        return new(ββ_1)
    end
end

function interaction_over_grid!(interaction, x1, grid, V::CalogeroSutherland)
    interaction .= V.ββ_1 ./ ((grid .- x1).^2 .+ 0.1)
end

abstract type NonInteracting <: Interaction end

struct HarmonicOscillator <: NonInteracting
    ω2::Float64
    HarmonicOscillator(ω) = new(ω^2)
end