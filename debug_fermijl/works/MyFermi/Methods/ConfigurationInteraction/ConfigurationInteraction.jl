"""
    MyFermi.ConfigurationInteraction

Module for running ConfigurationInteraction computations in MyFermi.
"""
module ConfigurationInteraction
# Import MyFermi basics
using MyFermi
using MyFermi.Options
using MyFermi.Integrals
using MyFermi.Orbitals

function ci_header()
    banner = 
raw"""
================================================================================
//                       Configuration Interaction                            \\     
//                        Module by G.J.R. Aroeira                            \\       
================================================================================
"""
    output(banner)
end

"""
    MyFermi.ConfigurationInteraction.AbstractCIWavefunction

MyFermi abstract type common to all Configuration Interaction wavefunctions

_struct tree:_

**AbstractCIWavefunction** <: AbstractWavefunction
"""
abstract type AbstractCIWavefunction <: MyFermi.AbstractWavefunction end

include("FCI/FCI.jl")

end #module CI
