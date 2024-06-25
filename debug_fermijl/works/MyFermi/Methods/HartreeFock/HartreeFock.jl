"""
    MyFermi.HartreeFock

Module for running Hartree--Fock computations in MyFermi.

# Methods

    > MyFermi.HartreeFock.RHF
    > MyFermi.HartreeFock.UHF
"""
module HartreeFock
# Import MyFermi basics
using LinearAlgebra: hermitian, vcat
using MyFermi
using MyFermi.Options
using MyFermi.Integrals
using MyFermi.Orbitals

function hf_header()
    output(repeat("=",80))
    output("|{:33}{:^12}{:33}|", "", "Hartree-Fock", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:25}{:^28}{:25}|", "", "G.J.R Aroeira and M.M. Davis", "")
    output(repeat("=",80))
end

function uhf_header()
    output(repeat("=",80))
    output("|{:33}{:^12}{:33}|", "", "Hartree-Fock", "")
    output("|{:34}{:^9}{:34}|", "", "Module  by","")
    output("|{:17}{:^44}{:17}|", "", "G.J.R Aroeira, M.M. Davis, and S.M. Goodlett", "")
    output(repeat("=",80))
end

"""
    MyFermi.HartreeFock.AbstractHFWavefunction

Abstract type common to all Hartree-Fock wave functions.

Struct tree

**AbstractHFWavefunction** <: AbstractWavefunction
"""
abstract type AbstractHFWavefunction <: MyFermi.AbstractWavefunction end

# Different Hartree-Fock methods are included here:
# Restricted Hartree--Fock
include("RHF/RHF.jl")

# UHF
include("UHF/UHF.jl")

end #module
