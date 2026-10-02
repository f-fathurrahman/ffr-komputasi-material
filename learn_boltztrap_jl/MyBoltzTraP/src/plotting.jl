# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Plotting stubs
#
# Implementations moved to ext/BoltzTraPPlotsExt.jl

"""
    plot_bands(result::InterpolationResult; kwargs...)
    plot_bands(file::String; kwargs...)

Plot band structure. Requires `using Plots`.
"""
function plot_bands end

"""
    plot_transport(result::TransportResult; kwargs...)
    plot_transport(file::String; kwargs...)

Plot transport properties. Requires `using Plots`.
"""
function plot_transport end

"""
    savefig_publication(bpd::BandPlotData, path; kwargs...)
    savefig_publication(bpds::AbstractVector{<:BandPlotData}, path; kwargs...)
    savefig_publication(tpd::TransportPlotData, path; kwargs...)
    savefig_publication(tpds::AbstractVector{<:TransportPlotData}, path; kwargs...)

Save a publication-quality figure via the PythonPlot extension.
Requires `using PythonPlot` to activate `BoltzTraPPythonPlotExt`.

The output format is inferred from the file extension (`.pdf`, `.png`, `.svg`, ...).

# Keyword arguments

  - `axis_width_cm = 8.0`, `axis_height_cm = 6.0` — Axis size in centimeters.
  - `ylims = nothing` — Y-axis range (auto if `nothing`; for `BandPlotData` defaults to `(emin, emax)`).
  - `title = ""` — Plot title (suppressed if empty; falls back to `data.title`).
  - `layout = (1, length(data))` — Subplot grid for the vector overloads.
  - `overlay = false` — Plot all `bpds` on a single axis (vector `BandPlotData` overload only).

See also: [`plot_bands`](@ref), [`plot_transport`](@ref).
"""
function savefig_publication end
