# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Plotting data types

"""
    BandPlotData

Precomputed band structure data for plotting, independent of plotting library.

# Fields

  - `kdist::Vector{Float64}` — cumulative k-path distances
  - `ebands_eV::Matrix{Float64}` — band energies in eV relative to Fermi (nbands × nk)
  - `unique_labels::Vector{String}` — high-symmetry point labels (merged with "|")
  - `unique_positions::Vector{Float64}` — high-symmetry point positions
  - `emin::Float64` — display energy minimum (eV)
  - `emax::Float64` — display energy maximum (eV)
  - `fermi_Ha::Float64` — Fermi energy in Hartree
  - `title::String` — plot title
"""
struct BandPlotData
    kdist::Vector{Float64}
    ebands_eV::Matrix{Float64}
    unique_labels::Vector{String}
    unique_positions::Vector{Float64}
    emin::Float64
    emax::Float64
    fermi_Ha::Float64
    title::String
end

"""
    TransportPlotData

Precomputed transport property data for plotting, independent of plotting library.

# Fields

  - `xdata::Vector{Float64}` — x-axis data (mu-Ef in eV, or temperature in K)
  - `ydata::Vector{Float64}` — y-axis data (already unit-converted)
  - `xlabel::String` — x-axis label
  - `ylabel::String` — y-axis label
  - `title::String` — plot title
"""
struct TransportPlotData
    xdata::Vector{Float64}
    ydata::Vector{Float64}
    xlabel::String
    ylabel::String
    title::String
end

# ── TransportPlotData builder ──
#
# Centralizes label formatting, unit conversion, and input validation so all
# downstream backends (RecipesBase / Plots / PythonPlot) consume the same
# normalized data. Pure: depends only on LinearAlgebra and HA_TO_EV.

# Component string → tensor index
const _COMPONENT_INDEX = Dict(
    "xx" => (1, 1),
    "yy" => (2, 2),
    "zz" => (3, 3),
    "xy" => (1, 2),
    "xz" => (1, 3),
    "yz" => (2, 3),
    "yx" => (2, 1),
    "zx" => (3, 1),
    "zy" => (3, 2),
)

#=
    _QuantitySpec(symbol, unit, factor, field)

Display metadata for a transport quantity. Internal.
=#
struct _QuantitySpec
    symbol::String   # display, e.g. "S", "σ", "κ"
    unit::String     # e.g. "μV/K"
    factor::Float64  # multiplicative conversion (V/K → μV/K = 1e6)
    field::Symbol    # TransportResult field name
end

const _QUANTITY_SPECS = Dict(
    "seebeck" => _QuantitySpec("S", "μV/K", 1e6, :seebeck),
    "S" => _QuantitySpec("S", "μV/K", 1e6, :seebeck),
    "sigma" => _QuantitySpec("σ", "S/m", 1.0, :sigma),
    "conductivity" => _QuantitySpec("σ", "S/m", 1.0, :sigma),
    "kappa" => _QuantitySpec("κ", "W/m/K", 1.0, :kappa),
    "thermal" => _QuantitySpec("κ", "W/m/K", 1.0, :kappa),
)

"""
    build_transport_plot_data(result::TransportResult;
        quantity="seebeck", component="xx", abscissa="mu",
        temperature=300.0, mu_index=0) -> TransportPlotData

Build a `TransportPlotData` from a `TransportResult` for plotting via any
backend (RecipesBase recipe, Plots, PythonPlot). Performs label formatting,
unit conversion, and input validation in a single place.

# Keyword arguments

  - `quantity::AbstractString = "seebeck"` — `"seebeck"`/`"S"`, `"sigma"`/`"conductivity"`, or `"kappa"`/`"thermal"`.
  - `component::AbstractString = "xx"` — Tensor component (`xx`, `yy`, `zz`, `xy`, `xz`, `yz`, `yx`, `zx`, `zy`).
  - `abscissa::AbstractString = "mu"` — `"mu"` for vs chemical potential, `"T"` for vs temperature.
  - `temperature::Float64 = 300.0` — Temperature [K] used when `abscissa = "mu"`.
  - `mu_index::Int = 0` — μ-grid index used when `abscissa = "T"`. `0` selects the value closest to the Fermi level.

Throws `ArgumentError` for unknown `quantity`, `component`, or `abscissa`,
or when `temperature` is not present in `result.temperatures`.
"""
function build_transport_plot_data(
    result::TransportResult;
    quantity::AbstractString = "seebeck",
    component::AbstractString = "xx",
    abscissa::AbstractString = "mu",
    temperature::Float64 = 300.0,
    mu_index::Int = 0,
)
    qspec = get(_QUANTITY_SPECS, quantity, nothing)
    isnothing(qspec) && throw(
        ArgumentError(
            "Unknown quantity=$(repr(quantity)). Expected one of $(sort!(collect(keys(_QUANTITY_SPECS)))).",
        ),
    )
    comp_lc = lowercase(component)
    haskey(_COMPONENT_INDEX, comp_lc) || throw(
        ArgumentError(
            "Unknown component=$(repr(component)). Expected xx/yy/zz/xy/xz/yz/yx/zx/zy.",
        ),
    )
    abscissa in ("mu", "T") || throw(
        ArgumentError("Unknown abscissa=$(repr(abscissa)). Expected \"mu\" or \"T\"."),
    )

    i, j = _COMPONENT_INDEX[comp_lc]
    data = getfield(result, qspec.field)  # (3, 3, nT, nμ)

    md = result.metadata
    fermi_Ha =
        haskey(md, "fermi_dft_Ha") ? md["fermi_dft_Ha"] :
        haskey(md, "fermi_Ha") ? md["fermi_Ha"] : 0.0
    fermi_eV = fermi_Ha * HA_TO_EV

    if abscissa == "mu"
        iT = findfirst(t -> isapprox(t, temperature, atol = 1.0), result.temperatures)
        isnothing(iT) && throw(
            ArgumentError(
                "Temperature $(temperature) K not found in result.temperatures = $(result.temperatures).",
            ),
        )
        ydata = data[i, j, iT, :] .* qspec.factor
        xdata = result.mu_values .- fermi_eV
        xlabel = "μ - E_F (eV)"
        title = "T = $(Int(temperature)) K"
    else  # "T"
        if mu_index == 0
            mu_index = argmin(abs.(result.mu_values .- fermi_eV))
        end
        ydata = data[i, j, :, mu_index] .* qspec.factor
        xdata = result.temperatures
        xlabel = "Temperature (K)"
        title = "μ = $(round(result.mu_values[mu_index] - fermi_eV, digits = 3)) eV"
    end

    ylabel = "$(qspec.symbol)_$(comp_lc) ($(qspec.unit))"

    return TransportPlotData(xdata, ydata, xlabel, ylabel, title)
end
