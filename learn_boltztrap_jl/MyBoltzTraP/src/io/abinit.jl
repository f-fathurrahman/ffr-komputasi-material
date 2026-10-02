# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
#=
ABINIT NetCDF file loader.

This module provides a stub for the ABINIT loader. The actual implementation
is in the NCDatasets extension (ext/BoltzTraPNCDatasetsExt.jl).

ABINIT outputs data in NetCDF format (*_GSR.nc files), which requires
NCDatasets.jl to read.
=#

export load_abinit

"""
    load_abinit(directory::String) -> [`DFTData`](@ref)

Load DFT data from ABINIT NetCDF output.

# Arguments

  - `directory`: Path to directory containing `*_GSR.nc` file

# Returns

  - [`DFTData`](@ref): DFT calculation data

# Notes

  - Requires NCDatasets.jl extension
  - Install with: `using Pkg; Pkg.add("NCDatasets")`

# File Format

ABINIT outputs Ground State Results in NetCDF format:

  - `*_GSR.nc`: Contains eigenvalues, k-points, structure, Fermi energy
"""
function load_abinit(args...)
    error(
        "load_abinit requires NCDatasets.jl. Run `using NCDatasets` before `using BoltzTraP`.",
    )
end

#=
    _detect_abinit(directory::String) -> Bool

Detect if directory contains ABINIT output files.

Returns true if a *_GSR.nc file is found.
=#
function _detect_abinit(directory::String)
    !isdir(directory) && return false
    for f in readdir(directory)
        if endswith(f, "_GSR.nc") && isfile(joinpath(directory, f))
            return true
        end
    end
    return false
end
