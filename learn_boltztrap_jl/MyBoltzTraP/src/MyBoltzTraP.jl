# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl

"""
    BoltzTraP

Julia implementation of BoltzTraP2, a band structure interpolation
and semi-classical transport coefficient calculator.

Based on BoltzTraP2 by G. K. H. Madsen, J. Carrete, and M. J. Verstraete.
"""
module MyBoltzTraP

using LinearAlgebra
using StaticArrays
using FFTW
using StatsBase
using Spglib

# Units and constants
include("units.jl")

# Spin type definitions (must be before types.jl for BandData{ST})
include("spintypes.jl")

# Type definitions
include("types.jl")

# Sphere module (bounds, rotations, tensor basis)
include("sphere.jl")

# Symmetry operations (Spglib wrapper)
include("symmetry.jl")

# Equivalence classes for reciprocal lattice points
include("equivalences.jl")

# Fourier interpolation
include("interpolation.jl")

# BZ sampling abstractions (dispatch axis for solve_transport)
include("sampling.jl")

# Band reconstruction via FFT
include("reconstruction.jl")

# Band library (DOS, Fermi integrals, Onsager coefficients)
include("bandlib.jl")

# I/O modules
include("io/serialization.jl")
include("io/vasp.jl")
include("io/qe.jl")
include("io/wien2k.jl")
include("io/gene.jl")
include("io/abinit.jl")
include("io/loader.jl")  # Auto-detection (must be after format-specific loaders)

# High-level workflow
include("workflow.jl")

# Method-dispatched transport API (AbstractBZSampling axis)
include("solve_transport.jl")

# Plotting types (backend-independent)
include("plottypes.jl")

# Plotting functions
include("plotting.jl")
include("fallbacks.jl")

# Describe functions (result inspection)
include("describe.jl")

# Extension stubs (implemented in ext/ when dependencies are loaded)
"""
    load_dftk(scfres) -> DFTData{1}

Extract band structure data from DFTK SCF result.

This function requires DFTK.jl to be loaded. Load it with:

```julia
using DFTK
using BoltzTraP
```

See the DFTK extension documentation for details.
"""
function load_dftk(args...)
    error("load_dftk requires DFTK.jl. Run `using DFTK` before `using BoltzTraP`.")
end

# CLI (Comonicon.jl)
include("cli.jl")

# =============================================================================
# Public API (Minimal)
# =============================================================================

# Workflow
export run_interpolate, run_integrate
export solve_transport

# Types
export InterpolationResult, TransportResult
export DFTData, NonMagneticData, SpinPolarizedData, nspin, is_magnetic

# BZ sampling + interpolators (the two dispatch axes of solve_transport)
export AbstractBZSampling, UniformMesh, TransportSystem
export AbstractInterpolator, FourierInterpolator, WannierInterpolator
export interpolate_velocities

# Spin types (v0.3 magnetic support)
export SpinType, Unpolarized, Collinear, NonCollinear
export BandData

# I/O
export load_dft, load_vasp, load_qe, load_wien2k, load_gene, load_abinit, load_dftk
export detected_format
export save_interpolation, load_interpolation
export save_integrate, load_integrate

# Plotting
export plot_bands, plot_transport
export BandPlotData, TransportPlotData
export build_transport_plot_data
export savefig_publication

# Describe (result inspection)
export describe

# Band library (additional transport functions)
export calc_cv
export apply_scissor

# =============================================================================
# Internal Functions
# =============================================================================
# The following functions are not part of the public API and are subject to
# change without notice. They are used internally by the workflow functions
# but are not intended for direct use.
#
# - `calc_N`: Computes the electron count for a given chemical potential.
# - `solve_for_mu`: Solves for the chemical potential for a given electron count.
# =============================================================================

end # module BoltzTraP
