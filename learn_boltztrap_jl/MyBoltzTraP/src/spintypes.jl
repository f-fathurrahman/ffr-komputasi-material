# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - SpinType definitions for magnetic material support

"""
    SpinType

Abstract type representing the spin polarization of a material.

Subtypes:

  - [`Unpolarized`](@ref): Non-spin-polarized (nonmagnetic)
  - [`Collinear`](@ref): Collinear spin-polarized (magnetic, scalar magnetization)
  - [`NonCollinear`](@ref): Non-collinear spin-polarized (magnetic, vector magnetization)
"""
abstract type SpinType end

"""
    Unpolarized <: SpinType

Spin type for non-spin-polarized (nonmagnetic) calculations.

Properties:

  - `dosweight = 2.0` (spin degeneracy factor)
  - `nspinor = 1` (single spin channel)
"""
struct Unpolarized <: SpinType end

"""
    Collinear <: SpinType

Spin type for collinear spin-polarized calculations.

Properties:

  - `dosweight = 1.0` (no spin degeneracy)
  - `nspinor = 2` (two spin channels: up and down)
  - Magnetization is a scalar per atom
"""
struct Collinear <: SpinType end

"""
    NonCollinear <: SpinType

Spin type for non-collinear spin-polarized calculations.

Properties:

  - `dosweight = 1.0` (no spin degeneracy)
  - `nspinor = 1` (spinor treated together)
  - Magnetization is a 3D vector per atom
"""
struct NonCollinear <: SpinType end

# ============================================================================
# Utility functions
# ============================================================================

"""
    dosweight(::Type{<:SpinType}) -> Float64

Return the DOS weight (spin degeneracy factor) for the given spin type.

  - `Unpolarized`: 2.0 (both spins contribute equally)
  - `Collinear`: 1.0 (spins counted separately)
  - `NonCollinear`: 1.0 (spins counted separately)
"""
dosweight(::Type{Unpolarized}) = 2.0
dosweight(::Type{Collinear}) = 1.0
dosweight(::Type{NonCollinear}) = 1.0

"""
    nspinor(::Type{<:SpinType}) -> Int

Return the number of spin channels for the given spin type.

  - `Unpolarized`: 1 (single channel, degeneracy handled via dosweight)
  - `Collinear`: 2 (spin-up and spin-down channels)
  - `NonCollinear`: 1 (spinor components treated together)
"""
nspinor(::Type{Unpolarized}) = 1
nspinor(::Type{Collinear}) = 2
nspinor(::Type{NonCollinear}) = 1
