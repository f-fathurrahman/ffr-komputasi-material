abstract type Basis end
abstract type SpatialBasis <: Basis end

struct HarmonicOscillatorBasis <: SpatialBasis
    l::Int64  # number of basis functions, degrees?
    ω::Float64  # strength of harmonic oscillator potential
    hermites::Vector{Float64}
    hos::Vector{Float64}
    ho_der::Vector{Float64}
    ho_dder::Vector{Float64}
end

# Represent basis for one spatial point?
function HarmonicOscillatorBasis(l, ω)
    hermites = zeros(Float64, l)
    return HarmonicOscillatorBasis(
        l, ω, hermites,
        zeros(Float64, l), zeros(Float64, l), zeros(Float64, l)
    )
end

# FIXME: use SpecialFunctions?

# using recursion relation
function evaluate_basis!(hos, x::Float64, ho::HarmonicOscillatorBasis)
    l = ho.l
    ω = ho.ω
    hermites = ho.hermites
    
    x = √ω * x # scale

    hermites[1] = 1.0
    ho_fac = (ω / π)^0.25 * exp(-x^2 / 2)
    hos[1] = ho_fac * hermites[1]
    #
    hermites[2] = 2*x
    ho_fac *= 1 / √2
    hos[2] = ho_fac * hermites[2]
    #
    for i in 3:l
        hermites[i] = 2x * hermites[i-1] - 2*(i - 2) * hermites[i-2]
        ho_fac *= 1 / sqrt( 2*(i - 1) )
        hos[i] = ho_fac * hermites[i]
    end
    return
end


function evaluate_on_grid(ho::HarmonicOscillatorBasis, xgrid::Vector{Float64})
    Nbasis = ho.l
    hos = ho.hos
    Npoints = size(xgrid, 1)
    spfs = [zeros(Float64, Npoints) for _ in 1:Nbasis]
    for ip in 1:Npoints
        evaluate_basis!(hos, xgrid[ip], ho)
        for j in 1:Nbasis
            spfs[j][ip] = hos[j]
        end
    end
    return spfs
end



# Rename to SpinOrbitalBasis ???
struct SpinBasis{T} <: SpatialBasis
    base::T # The basis with no spin
    l::Int64 # number of basis (2 times larger than spatial basis)
end

function SpinBasis(base::SpatialBasis)
    @info "Constructing SpinBasis"
    return SpinBasis{typeof(base)}(base, 2*base.l)
end

function evaluate_on_grid(basis::SpinBasis, xgrid)
    Npoints = length(xgrid)
    Nbasis = basis.l # including spin
    NbasisSpatial = floor(Int64, Nbasis/2)
    res = Vector{Vector{Float64}}(undef, Nbasis)
    for i in 1:Nbasis
        res[i] = zeros(Float64, Npoints)
    end
    nospin = zeros(Float64, NbasisSpatial) # temporary array
    # loop over spatial points
    for i in 1:Npoints
        evaluate_basis!(nospin, xgrid[i], basis.base) # the basis functions evaluated at x
        for j in 1:NbasisSpatial
            res[2*j-1][i] = nospin[j]
            res[2*j][i] = nospin[j]
        end
    end 
    return res
end