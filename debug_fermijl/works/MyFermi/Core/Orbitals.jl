module Orbitals

using MyFermi
using MyGaussianBasis

export AbstractOrbitals, AtomicOrbitals, AbstractRestrictedOrbitals, AbstractUnrestrictedOrbitals
export GeneralRestrictedOrbitals
export RHFOrbitals
export UHFOrbitals

"""
    MyFermi.AbstractOrbitals

Abstract type common to all orbitals

_struct tree:_

**AbstractOrbitals**  (Top level)
"""
abstract type AbstractOrbitals end

struct AtomicOrbitals <: AbstractOrbitals 
    basisset::BasisSet
end

function AtomicOrbitals(mol:: MyMolecule, basis::String)
    return AtomicOrbitals(BasisSet(basis, mol.atoms))
end

abstract type AbstractRestrictedOrbitals <: AbstractOrbitals end
abstract type AbstractUnrestrictedOrbitals <: AbstractOrbitals end

"""
    MyFermi.GeneralRestrictedOrbitals

Dummy orbital struct that can be used to hold any type of Restricted orbitals.
Can be used to read in custom orbitals as

    Orbs = GeneralRestrictedOrbitals(X)

where X is an Array object. The molecule and basis will be deduced from the current options.
Alternatively, one can pass these informations explicitly:

    Orbs = GeneralRestrictedOrbitals(X, molecule=mol, name="myorbitals", basis="cc-pvdz")

`name` and `basis` are Strings, whereas `mol` is a `MyFermi.MyMolecule` object.

# Fields

    name      String with a label for the orbital
    basis     String indicating the basis set used to construct the orbitals
    molecule  MyMolecule object for which the orbitals were constructed
    C         NxN AbstractArray with the AO(lines) → MO(orbitals) coefficients

_struct tree:_

**GeneralRestrictedOrbitals** <: AbstractOrbitals
"""
struct GeneralRestrictedOrbitals{T} <: AbstractRestrictedOrbitals 
    molecule::MyMolecule
    basis::String
    sd_energy::T
    C::AbstractArray{T,2}
end

function GeneralRestrictedOrbitals(C::AbstractArray{T,2}; mol=nothing, basis="undef", sd_energy=zero(T)) where T <: AbstractFloat

    mol === nothing ? mol = MyFermi.MyMolecule() : nothing
    basis == "undef" ? basis = MyFermi.Options.get("basis") : nothing

    GeneralRestrictedOrbitals{T}(mol, basis, sd_energy, C)
end

"""
    MyFermi.HartreeFock.RHFOrbitals

Struct holding information about Restricted Hartree--Fock orbitals

# Fields

    molecule   MyMolecule object associated with the orbitals
    basis      Basis set used to compute the orbitals
    eps        Orbital energies, i.e. diagonal of the Fock matrix
    C          Coefficients of the AO->MO transformation matrix
"""
struct RHFOrbitals <: AbstractRestrictedOrbitals
    molecule::MyMolecule
    basis::String
    eps::AbstractArray{Float64,1}
    sd_energy::Float64
    C::AbstractArray{Float64,2}
end

struct UHFOrbitals <: AbstractUnrestrictedOrbitals
    molecule::MyMolecule
    basis::String
    epsα::AbstractArray{Float64,1}
    epsβ::AbstractArray{Float64,1}
    sd_energy::Float64
    Cα::AbstractArray{Float64,2}
    Cβ::AbstractArray{Float64,2}
end


end # Module