

AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(
    system,
    calculator::MyACEmd.MyACEpotential;
    kwargs...
)
    return ace_energy(calculator, system; kwargs...)
end



AtomsCalculators.@generate_interface function AtomsCalculators.forces(
    system,
    calculator::MyACEmd.MyACEpotential;
    kwargs...
)
    return ace_forces(calculator, system; kwargs...)
end



AtomsCalculators.@generate_interface function AtomsCalculators.virial(
    system,
    calculator::MyACEmd.MyACEpotential;
    kwargs...
)
    return ace_virial(calculator, system; kwargs...)
end


function AtomsCalculators.energy_forces(system, calculator::MyACEpotential; kwargs...)
    return ace_energy_forces(calculator, system;  kwargs...)
end


function AtomsCalculators.energy_forces_virial(system, calculator::MyACEpotential; kwargs...)
    return ace_energy_forces_virial(calculator, system;  kwargs...)
end