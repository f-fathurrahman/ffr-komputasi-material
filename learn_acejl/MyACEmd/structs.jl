
const default_length = u"Å"
const default_energy = u"eV"

struct MyACEpotential{TE,TL,TC}
    potentials::Vector{MyACE1.AbstractCalculator}
    energy_unit::TE
    length_unit::TL
    cutoff_unit::TC
    function MyACEpotential(potentials; energy_unit=default_energy, length_unit=default_length, cutoff_unit=length_unit)
        @assert dimension(energy_unit) == dimension(u"J")
        @assert dimension(length_unit) == dimension(u"m")
        @assert dimension(cutoff_unit) == dimension(u"m")
        new{typeof(energy_unit), typeof(length_unit), typeof(cutoff_unit)}(potentials, energy_unit, length_unit, cutoff_unit)
    end
end


function MyACEpotential(
    basis::MyACE1.MLIPs.IPBasis,
    constants, 
    oneparticle::Union{MyACE1.OneBody, Nothing}=nothing; 
    kwargs...
)
    if isnothing(oneparticle)
        pot = MyACE1.MLIPs.combine(basis, constants)
    else
        pot = MyACE1.MLIPs.SumIP(oneparticle, MyACE1.MLIPs.combine(basis, constants))
    end
    return MyACEpotential(pot.components; kwargs...)
end


function Base.iterate(acep::MyACEpotential, state::Int=1)
    if 0 < state <= length( acep )
        return acep.potentials[state], state + 1
    else
        return nothing
    end
end

Base.length(acep::MyACEpotential) = length(acep.potentials)

Base.getindex(acep::MyACEpotential, i) = acep.potentials[i]

function Base.show(io::IO, ::MIME"text/plain", acep::MyACEpotential)
    print(io, "MyACE potential consisting of ", length(acep), " subpotentials")
end