# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""
    AbstractBZSampling

Abstract type for Brillouin-zone sampling strategies used by transport
calculations. Concrete subtypes (e.g. `UniformMesh`) carry the parameters
that determine how `solve_transport` integrates over the BZ.

This is the dispatch axis for BZ sampling strategies; alternative methods
(adaptive grids, AutoBZ-style integration, wavelet meshes) can be added
later as additional concrete subtypes.
"""
abstract type AbstractBZSampling end

"""
    UniformMesh(nk=nothing)
    UniformMesh((nk1, nk2, nk3))

Uniform Monkhorst-Pack-style mesh. The `nk` field asks for a mesh density:

  - `nk === nothing`: let the interpolator choose the mesh size.
  - `nk::Tuple{Int,Int,Int}`: an explicit mesh size.

Whether an explicit size is honoured depends on the interpolator. The Wannier
path requires one and raises an error without it. The Fourier path takes its
grid from the interpolation basis, where the star functions fix the smallest
grid that represents them exactly, and warns that an explicit size is ignored.

The grid actually used is recorded in `TransportResult.metadata["sampling_nk"]`.
"""
struct UniformMesh <: AbstractBZSampling
    nk::Union{Nothing,Tuple{Int,Int,Int}}
end

UniformMesh() = UniformMesh(nothing)

"""
    TransportSystem(lattice, fermi, nelect, dosweight, spintype)
    TransportSystem(interp::InterpolationResult)

System-level constants needed for transport calculations: lattice, Fermi
energy, electron count, DOS weight, and spin type. These quantities are
conventionally spread across `InterpolationResult.metadata` and the
underlying band data; this struct unifies them as a single, dispatched
data carrier for the new `solve_transport` API.

Fields:

  - `lattice::SMatrix{3,3,Float64,9}`: 3×3 lattice vectors (columns) in Bohr
  - `fermi::Float64`: Fermi energy in Hartree
  - `nelect::Float64`: number of electrons per unit cell
  - `dosweight::Float64`: DOS weight (2.0 = non-magnetic, 1.0 = spin-polarized)
  - `spintype::SpinType`: `Unpolarized()` / `Collinear()` / `NonCollinear()`

The `TransportSystem(::InterpolationResult)` constructor is defined in
`io/serialization.jl` (after `InterpolationResult` itself is defined).
"""
struct TransportSystem
    lattice::SMatrix{3,3,Float64,9}
    fermi::Float64
    nelect::Float64
    dosweight::Float64
    spintype::SpinType
end

"""
    WannierInterpolator(prefix::String; spintype=Unpolarized())

Wannier-based interpolator constructed from a Wannier90 file prefix
(`.chk` and related files). Concrete methods for `interpolate_bands`,
`interpolate_velocities`, and the `solve_transport` dispatch are provided
by `BoltzTraPWannierExt`, which loads when both `Wannier` and `BoltzTraP`
are imported.

The `model` field is held as `Any` because `Wannier.InterpModel` is
defined in the optional `Wannier.jl` package; the stub constructor below
errors when `using Wannier` has not been invoked, ensuring users see a
clear message rather than a `MethodError`.

Fields:

  - `model::Any`: a `Wannier.InterpModel` instance at runtime
  - `lattice::SMatrix{3,3,Float64,9}`: 3×3 lattice vectors (columns) in
    Bohr (converted from Wannier's Ångström convention by the extension)
  - `spintype::ST`: `Unpolarized()` / `Collinear()` / `NonCollinear()`
"""
struct WannierInterpolator{ST<:SpinType} <: AbstractInterpolator
    model::Any
    lattice::SMatrix{3,3,Float64,9}
    spintype::ST
end

"""
    TransportSystem(wi::WannierInterpolator; fermi, nelect, spintype, dosweight)

Build a [`TransportSystem`](@ref) from a `WannierInterpolator`.

The lattice and the spin type are taken from `wi`. The Fermi energy and the
electron count are not present in the Wannier90 files, so they have no default
and must be supplied by the caller.

# Keyword Arguments

  - `fermi`: Fermi energy in Ha (required)
  - `nelect`: number of electrons per unit cell (required)
  - `spintype=wi.spintype`: spin type carried by the interpolator
  - `dosweight`: DOS weight; defaults to the spin degeneracy implied by
    `spintype`, i.e. `2.0` for `Unpolarized()` and `1.0` otherwise
"""
function TransportSystem(
    wi::WannierInterpolator;
    fermi::Real,
    nelect::Real,
    spintype::SpinType = wi.spintype,
    dosweight::Real = (spintype isa Unpolarized ? 2.0 : 1.0),
)
    return TransportSystem(wi.lattice, fermi, nelect, dosweight, spintype)
end

# Stub: file-load constructor; the real method lives in BoltzTraPWannierExt.
# The signature is `(args...; kwargs...)` — less specific than the extension's
# `(prefix::String; spintype)` — so Julia dispatches to the extension method
# when Wannier.jl is loaded and falls back to this error otherwise. This
# matches the load_dftk / load_abinit weakdep pattern used elsewhere in BT.jl.
function WannierInterpolator(args...; kwargs...)
    error(
        "WannierInterpolator(::String) requires Wannier.jl to be loaded. " *
        "Run `using Wannier` before constructing a WannierInterpolator from " *
        "a Wannier90 prefix.",
    )
end
