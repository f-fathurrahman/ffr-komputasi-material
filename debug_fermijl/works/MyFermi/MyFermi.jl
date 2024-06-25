"""
    MyFermi Quantum Chemistry Module

Authors: G. J. R. Aroeira, M. M. Davis, J. M. Turney, and H. F. Schaefer

GitHub: [MyFermi.jl](https://github.com/MyFermiQC/MyFermi.jl)
"""
module MyFermi

using LinearAlgebra: Threads
"""
    MyFermi.AbstractWavefunction

Abstract type common to all wave functions.

_struct tree:_

**AbstractWavefunction**  (Top level)
"""
abstract type AbstractWavefunction end

include("Core/Options.jl")                             
include("Backend/Arrays.jl")
include("Core/Output.jl")
include("Backend/PhysicalConstants.jl")
include("Core/DIIS.jl")
include("Core/Molecule.jl")                               
include("Core/Orbitals.jl")
include("Core/Integrals/IntegralHelper.jl")
include("Methods/HartreeFock/HartreeFock.jl")
include("Methods/MollerPlesset/MollerPlesset.jl")
include("Methods/CoupledCluster/CoupledCluster.jl")
include("Methods/ConfigurationInteraction/ConfigurationInteraction.jl")

include("Tools/SinglePointEnergy.jl")
include("Tools/FiniteDifferences.jl")
include("Tools/Gradients.jl")

end # module
