# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using JLD2
using JSON3
using CodecXz
using Dates

# ============================================================================
# Data Structures for Serialization
# ============================================================================

"""
    InterpolationResult

Container for interpolation results from [`run_interpolate`](@ref).

# Fields

  - `coeffs`: Fourier coefficients (nbands × neq)
  - `equivalences`: Vector of equivalence class matrices
  - `lattvec`: 3×3 lattice vectors (columns, in Bohr)
  - `atoms`: Atomic structure (positions, types) - optional
  - `metadata`: Dictionary with additional information

# Metadata Keys

  - `fermi`: Fermi energy in Ha
  - `nelect`: Number of electrons
  - `dosweight`: DOS weight (2.0 for non-spin-polarized, 1.0 for collinear)
  - `spintype`: Spin type ("Unpolarized" or "Collinear")
  - `selected_bands`: Band indices used (UnitRange)
  - `spacegroup_number`: International space group number (1-230)
  - `spacegroup_symbol`: Space group symbol (e.g., "Fd-3m")
  - `source`: Original DFT file path
  - `created`: ISO 8601 timestamp
  - `generator`: "BoltzTraP.jl"

See also: [`run_interpolate`](@ref), [`run_integrate`](@ref), [`save_interpolation`](@ref), [`load_interpolation`](@ref)
"""
struct InterpolationResult
    coeffs::Matrix{ComplexF64}
    equivalences::Vector{Matrix{Int}}
    lattvec::Matrix{Float64}
    atoms::Union{Nothing,Dict{String,Any}}
    metadata::Dict{String,Any}
end

"""
    InterpolationResult(interp::FourierInterpolator; atoms=nothing, metadata=Dict())

Create InterpolationResult from a FourierInterpolator.
"""
function InterpolationResult(
    interp::FourierInterpolator;
    atoms = nothing,
    metadata = Dict{String,Any}(),
)
    # Convert equivalences to plain Matrix{Int} for serialization
    equivs = [Matrix{Int}(eq) for eq in interp.equivalences]
    lattvec = Matrix{Float64}(interp.lattvec)

    meta = Dict{String,Any}(
        "version" => "0.1.0",
        "created" => string(now()),
        "nbands" => size(interp.coeffs, 1),
        "neq" => length(interp.equivalences),
        "spintype" => "Unpolarized",  # v0.2 forward compatibility
    )
    merge!(meta, metadata)

    return InterpolationResult(interp.coeffs, equivs, lattvec, atoms, meta)
end

"""
    TransportSystem(interp::InterpolationResult) -> TransportSystem

Build a [`TransportSystem`](@ref) from an `InterpolationResult`, extracting
lattice (Bohr), Fermi energy (Ha), electron count, DOS weight, and spin type
from the standard metadata keys.

Required metadata keys: `"fermi"` (Ha), `"nelect"`, `"dosweight"`,
`"spintype"` (`"Unpolarized"` / `"Collinear"` / `"NonCollinear"`).
"""
function TransportSystem(interp::InterpolationResult)
    lat = SMatrix{3,3,Float64}(interp.lattvec)
    fermi = interp.metadata["fermi"]
    nelect = interp.metadata["nelect"]
    dosweight = interp.metadata["dosweight"]
    # v0.1 metadata (pre-spintype era) defaults to Unpolarized.
    spintype_str = get(interp.metadata, "spintype", "Unpolarized")
    spintype = if spintype_str == "Unpolarized"
        Unpolarized()
    elseif spintype_str == "Collinear"
        Collinear()
    elseif spintype_str == "NonCollinear"
        NonCollinear()
    else
        error("TransportSystem: unknown spintype \"$spintype_str\" in metadata")
    end
    return TransportSystem(lat, fermi, nelect, dosweight, spintype)
end

"""
    FourierInterpolator(interp::InterpolationResult) -> FourierInterpolator

Build a `FourierInterpolator` from an `InterpolationResult`, reusing the
stored Fourier coefficients, equivalence classes, and lattice vectors.

The inverse direction is [`InterpolationResult`](@ref)`(::FourierInterpolator)`,
which bundles an interpolator together with atoms and metadata for
serialization.
"""
function FourierInterpolator(interp::InterpolationResult)
    return FourierInterpolator(
        interp.coeffs,
        interp.equivalences,
        SMatrix{3,3,Float64}(interp.lattvec),
    )
end

"""
    TransportResult

Container for transport coefficients from [`run_integrate`](@ref).

Contains transport coefficients and other quantities computed by BZ integration.

# Fields

  - `temperatures`: Temperature array in Kelvin (K)
  - `mu_values`: Chemical potential array in eV
  - `sigma`: Electrical conductivity/τ tensor (3×3×nT×nμ)
  - `seebeck`: Seebeck coefficient tensor (3×3×nT×nμ) in V/K
  - `kappa`: Electronic thermal conductivity/τ tensor (3×3×nT×nμ)
  - `dos`: Density of states data (optional Dict with `epsilon`, `dos`, `vvdos`)
  - `metadata`: Dictionary with additional information

# Metadata Keys

Stamped by [`solve_transport`](@ref):

  - `interpolator`: `"Fourier"` or `"Wannier"`
  - `sampling`: name of the BZ sampling type, e.g. `"UniformMesh"`
  - `nelect`, `dosweight`, `spintype`: system constants used for the integration
  - `fermi_Ha`, `fermi_eV`: Fermi energy of the interpolated bands
  - `nbands`, `npts_fft`, `npts_dos`: grid sizes
  - further provenance keys, and `scissor_eV` / `scissor_Ha` when a scissor
    correction is applied

Carried over from the interpolation metadata by [`run_integrate`](@ref):

  - `source`: Original DFT file path

See also: [`run_integrate`](@ref), [`save_integrate`](@ref), [`load_integrate`](@ref)
"""
struct TransportResult
    temperatures::Vector{Float64}
    mu_values::Vector{Float64}
    sigma::Array{Float64,4}
    seebeck::Array{Float64,4}
    kappa::Array{Float64,4}
    dos::Union{Nothing,Dict{String,Any}}
    metadata::Dict{String,Any}
end

# ============================================================================
# High-level API
# ============================================================================

"""
    save_interpolation(filename, result::InterpolationResult)
    save_interpolation(filename, interp::FourierInterpolator; kwargs...)

Save interpolation results to file.

Supported formats:

  - `.jld2`: Julia native (HDF5-based, fast)
  - `.bt2`: Python BoltzTraP2 compatible (JSON + LZMA)
"""
function save_interpolation(filename::String, result::InterpolationResult)
    ext = splitext(filename)[2]

    if ext == ".jld2"
        save_interpolation_jld2(filename, result)
    elseif ext == ".bt2"
        save_interpolation_bt2(filename, result)
    else
        error("Unsupported file format: $ext. Use .jld2 or .bt2")
    end
end

function save_interpolation(filename::String, interp::FourierInterpolator; kwargs...)
    result = InterpolationResult(interp; kwargs...)
    save_interpolation(filename, result)
end

"""
    load_interpolation(filename) -> InterpolationResult

Load interpolation results from file.
"""
function load_interpolation(filename::String)
    ext = splitext(filename)[2]

    if ext == ".jld2"
        return load_interpolation_jld2(filename)
    elseif ext == ".bt2"
        return load_interpolation_bt2(filename)
    else
        error("Unsupported file format: $ext")
    end
end

"""
    save_integrate(filename, result::TransportResult)

Save integration results to file.

Supported formats:

  - `.jld2`: Julia native (full data)
  - `.csv`: Text format (scalar averages only)
"""
function save_integrate(filename::String, result::TransportResult)
    ext = splitext(filename)[2]

    if ext == ".jld2"
        save_integrate_jld2(filename, result)
    elseif ext == ".csv"
        save_integrate_csv(filename, result)
    else
        error("Unsupported file format: $ext. Use .jld2 or .csv")
    end
end

"""
    load_integrate(filename) -> TransportResult

Load integration results from file.
"""
function load_integrate(filename::String)
    ext = splitext(filename)[2]

    if ext == ".jld2"
        return load_integrate_jld2(filename)
    else
        error("Unsupported file format: $ext for loading. Use .jld2")
    end
end

# ============================================================================
# JLD2 Format (Julia Native)
# ============================================================================

function save_interpolation_jld2(filename::String, result::InterpolationResult)
    jldsave(
        filename;
        coeffs = result.coeffs,
        equivalences = result.equivalences,
        lattvec = result.lattvec,
        atoms = result.atoms,
        metadata = result.metadata,
    )
end

function load_interpolation_jld2(filename::String)
    data = load(filename)
    metadata = get(data, "metadata", Dict{String,Any}())

    # v0.1 backward compatibility: add spintype default if missing
    if !haskey(metadata, "spintype")
        metadata["spintype"] = "Unpolarized"
    end

    return InterpolationResult(
        data["coeffs"],
        data["equivalences"],
        data["lattvec"],
        get(data, "atoms", nothing),
        metadata,
    )
end

function save_integrate_jld2(filename::String, result::TransportResult)
    jldsave(
        filename;
        temperatures = result.temperatures,
        mu_values = result.mu_values,
        sigma = result.sigma,
        seebeck = result.seebeck,
        kappa = result.kappa,
        dos = result.dos,
        metadata = result.metadata,
    )
end

function load_integrate_jld2(filename::String)
    data = load(filename)
    metadata = get(data, "metadata", Dict{String,Any}())

    # v0.1 backward compatibility: add spintype default if missing
    if !haskey(metadata, "spintype")
        metadata["spintype"] = "Unpolarized"
    end

    return TransportResult(
        data["temperatures"],
        data["mu_values"],
        data["sigma"],
        data["seebeck"],
        data["kappa"],
        get(data, "dos", nothing),
        metadata,
    )
end

# ============================================================================
# BT2 Format (Python BoltzTraP2 Compatible)
# ============================================================================

#=
Convert InterpolationResult to Python-compatible dictionary format.
=#
function to_bt2_dict(result::InterpolationResult)
    # Convert complex coefficients to real/imag pairs for JSON
    coeffs_real = real.(result.coeffs)
    coeffs_imag = imag.(result.coeffs)

    Dict(
        "coeffs_real" => collect(eachrow(coeffs_real)),
        "coeffs_imag" => collect(eachrow(coeffs_imag)),
        "equivalences" => [collect(eachrow(eq)) for eq in result.equivalences],
        "lattvec" => collect(eachrow(result.lattvec)),
        "atoms" => result.atoms,
        "metadata" => result.metadata,
        "format_version" => "1.0",
        "generator" => "BoltzTraP.jl",
    )
end

#=
Convert Python-compatible dictionary to InterpolationResult.
=#
function from_bt2_dict(data::Dict)
    # Reconstruct complex coefficients
    coeffs_real = reduce(vcat, [reshape(row, 1, :) for row in data["coeffs_real"]])
    coeffs_imag = reduce(vcat, [reshape(row, 1, :) for row in data["coeffs_imag"]])
    coeffs = complex.(coeffs_real, coeffs_imag)

    # Reconstruct equivalences
    equivalences =
        [reduce(vcat, [reshape(row, 1, :) for row in eq]) for eq in data["equivalences"]]
    equivalences = [Matrix{Int}(eq) for eq in equivalences]

    # Reconstruct lattice vectors
    lattvec = reduce(vcat, [reshape(row, 1, :) for row in data["lattvec"]])
    lattvec = Matrix{Float64}(lattvec)

    atoms = get(data, "atoms", nothing)
    metadata = get(data, "metadata", Dict{String,Any}())

    # v0.1 backward compatibility: add spintype default if missing
    if !haskey(metadata, "spintype")
        metadata["spintype"] = "Unpolarized"
    end

    return InterpolationResult(coeffs, equivalences, lattvec, atoms, metadata)
end

function save_interpolation_bt2(filename::String, result::InterpolationResult)
    data = to_bt2_dict(result)
    json_str = JSON3.write(data)

    # Compress with XZ (LZMA2)
    compressed = transcode(XzCompressor, Vector{UInt8}(json_str))

    open(filename, "w") do io
        write(io, compressed)
    end
end

# ============================================================================
# Python BoltzTraP2 bt2 Format Support
# ============================================================================

#=
    _is_python_bt2_format(data) -> Bool

Check if data is in Python BoltzTraP2 format.
Python format: [data, equivalences, coeffs, metadata] array with BoltzTraP2_type markers.
=#
function _is_python_bt2_format(data)
    data isa AbstractVector &&
        length(data) == 4 &&
        data[1] isa AbstractDict &&
        get(data[1], "BoltzTraP2_type", "") == "DFTData"
end

#=
    _convert_python_bt2_array(obj) -> Array

Convert Python BoltzTraP2 Array type to Julia array.
Handles both Complex and regular arrays.
=#
function _convert_python_bt2_array(obj)
    if obj isa AbstractDict && get(obj, "BoltzTraP2_type", "") == "Array"
        if get(obj, "array_type", "") == "Complex"
            # Python stores as nested lists, need to convert to matrix
            real_data = obj["real"]
            imag_data = obj["imag"]
            # Convert nested vectors to matrix (nbands × neq)
            real_part =
                reduce(vcat, [reshape(collect(Float64, row), 1, :) for row in real_data])
            imag_part =
                reduce(vcat, [reshape(collect(Float64, row), 1, :) for row in imag_data])
            return complex.(real_part, imag_part)
        else
            # Regular array
            arr_data = obj["data"]
            return reduce(vcat, [reshape(collect(Float64, row), 1, :) for row in arr_data])
        end
    end
    return obj
end

#=
    _convert_python_equivalences(equiv_list) -> Vector{Matrix{Int}}

Convert Python equivalences list to Julia format.
Python stores each equivalence class as {"BoltzTraP2_type": "Array", "data": [[i,j,k], ...]}.
Optimized for performance with large equivalence lists.
=#
function _convert_python_equivalences(equiv_list)
    result = Vector{Matrix{Int}}()
    sizehint!(result, length(equiv_list))

    for eq in equiv_list
        bt2_type = get(eq, "BoltzTraP2_type", "")
        if bt2_type == "Array"
            eq_data = eq["data"]
            n = length(eq_data)
            if n > 0
                # Preallocate matrix for efficiency
                mat = Matrix{Int}(undef, n, 3)
                @inbounds for (j, point) in enumerate(eq_data)
                    mat[j, 1] = point[1]
                    mat[j, 2] = point[2]
                    mat[j, 3] = point[3]
                end
                push!(result, mat)
            else
                push!(result, Matrix{Int}(undef, 0, 3))
            end
        elseif eq isa AbstractVector
            # Fallback: direct array (shouldn't happen in Python bt2)
            n = length(eq)
            if n > 0
                mat = Matrix{Int}(undef, n, 3)
                @inbounds for (j, point) in enumerate(eq)
                    mat[j, 1] = point[1]
                    mat[j, 2] = point[2]
                    mat[j, 3] = point[3]
                end
                push!(result, mat)
            else
                push!(result, Matrix{Int}(undef, 0, 3))
            end
        else
            error("Unexpected equivalence format: $(typeof(eq))")
        end
    end
    return result
end

#=
    _json3_to_dict(obj) -> Dict{String,Any}

Convert JSON3 object to Dict{String,Any}, handling Symbol keys.
=#
function _json3_to_dict(obj)
    result = Dict{String,Any}()
    for (k, v) in pairs(obj)
        key = k isa Symbol ? String(k) : String(k)
        if v isa AbstractDict
            result[key] = _json3_to_dict(v)
        elseif v isa AbstractVector
            result[key] = [x isa AbstractDict ? _json3_to_dict(x) : x for x in v]
        else
            result[key] = v
        end
    end
    return result
end

#=
    _extract_cell_matrix(cell_data) -> Matrix{Float64}

Extract 3x3 cell matrix from JSON3 cell data.
Cell data may be wrapped in BoltzTraP2_type format.
=#
function _extract_cell_matrix(cell_data)
    # Check if cell is wrapped in BoltzTraP2_type format
    actual_data =
        if cell_data isa AbstractDict && get(cell_data, "BoltzTraP2_type", "") == "Array"
            cell_data["data"]
        else
            cell_data
        end

    # actual_data is a 3x3 nested array (each row is a lattice vector)
    rows = Vector{Vector{Float64}}()
    for row in actual_data
        push!(rows, Float64[x for x in row])
    end
    # Stack rows and transpose (Python row-major to Julia column-major)
    lattvec = reduce(vcat, [reshape(r, 1, :) for r in rows])
    return Matrix{Float64}(lattvec')
end

#=
    from_python_bt2(data::AbstractVector) -> InterpolationResult

Convert Python BoltzTraP2 format [data, equivalences, coeffs, metadata] to Julia.
=#
function from_python_bt2(data::AbstractVector)
    dft_data = data[1]
    equivalences = _convert_python_equivalences(data[2])
    coeffs = _convert_python_bt2_array(data[3])
    metadata = _json3_to_dict(data[4])

    # Extract lattice from ASE Atoms object
    atoms_dict = dft_data["atoms"]
    cell_data = atoms_dict["cell"]
    lattvec = _extract_cell_matrix(cell_data)

    # Convert atoms dict to Julia Dict
    atoms = _json3_to_dict(atoms_dict)

    # v0.1 backward compatibility: add spintype default if missing
    if !haskey(metadata, "spintype")
        metadata["spintype"] = "Unpolarized"
    end

    return InterpolationResult(coeffs, equivalences, lattvec, atoms, metadata)
end

function load_interpolation_bt2(filename::String)
    compressed = read(filename)
    json_bytes = transcode(XzDecompressor, compressed)
    json_str = String(json_bytes)

    # Parse JSON - use Any to handle both formats
    data = JSON3.read(json_str)

    # Auto-detect format
    if _is_python_bt2_format(data)
        return from_python_bt2(collect(data))
    else
        # Julia format - use _json3_to_dict for proper conversion
        return from_bt2_dict(_json3_to_dict(data))
    end
end

# ============================================================================
# CSV Format (Integration Results)
# ============================================================================

#=
    save_integrate_csv(filename, result::TransportResult)

Save integration results to CSV format.

Outputs scalar averages: σ_avg = (σ_xx + σ_yy + σ_zz) / 3
=#
function save_integrate_csv(filename::String, result::TransportResult)
    nT = length(result.temperatures)
    nμ = length(result.mu_values)

    open(filename, "w") do io
        # Header
        println(io, "# BoltzTraP.jl Transport Results")
        println(io, "# Generated: $(now())")
        if haskey(result.metadata, "material")
            println(io, "# Material: $(result.metadata["material"])")
        end
        println(io, "#")
        println(io, "# Columns:")
        println(io, "#   T: Temperature [K]")
        println(io, "#   mu: Chemical potential [eV]")
        println(io, "#   sigma: Electrical conductivity [S/m] (trace average)")
        println(io, "#   S: Seebeck coefficient [μV/K] (trace average)")
        println(io, "#   kappa_e: Electronic thermal conductivity [W/m/K] (trace average)")
        println(io, "#")
        println(io, "T,mu,sigma,S,kappa_e")

        for iT = 1:nT
            T = result.temperatures[iT]
            for iμ = 1:nμ
                μ = result.mu_values[iμ]

                # Trace averages (diagonal elements)
                σ_avg =
                    (
                        result.sigma[1, 1, iT, iμ] +
                        result.sigma[2, 2, iT, iμ] +
                        result.sigma[3, 3, iT, iμ]
                    ) / 3

                S_avg =
                    (
                        result.seebeck[1, 1, iT, iμ] +
                        result.seebeck[2, 2, iT, iμ] +
                        result.seebeck[3, 3, iT, iμ]
                    ) / 3
                S_avg_uV = S_avg * 1e6  # Convert V/K to μV/K

                κ_avg =
                    (
                        result.kappa[1, 1, iT, iμ] +
                        result.kappa[2, 2, iT, iμ] +
                        result.kappa[3, 3, iT, iμ]
                    ) / 3

                println(io, "$T,$μ,$σ_avg,$S_avg_uV,$κ_avg")
            end
        end
    end
end

#=
    save_integrate_trace_csv(filename, result::TransportResult)

Save only diagonal (trace) components of transport tensors to CSV.
More compact format for isotropic or cubic materials.
=#
function save_integrate_trace_csv(filename::String, result::TransportResult)
    nT = length(result.temperatures)
    nμ = length(result.mu_values)

    open(filename, "w") do io
        println(io, "# BoltzTraP.jl Transport Results (Diagonal Components)")
        println(io, "# Generated: $(now())")
        println(io, "#")
        println(
            io,
            "T,mu,sigma_xx,sigma_yy,sigma_zz,S_xx,S_yy,S_zz,kappa_xx,kappa_yy,kappa_zz",
        )

        for iT = 1:nT
            T = result.temperatures[iT]
            for iμ = 1:nμ
                μ = result.mu_values[iμ]

                σ_xx = result.sigma[1, 1, iT, iμ]
                σ_yy = result.sigma[2, 2, iT, iμ]
                σ_zz = result.sigma[3, 3, iT, iμ]

                S_xx = result.seebeck[1, 1, iT, iμ] * 1e6
                S_yy = result.seebeck[2, 2, iT, iμ] * 1e6
                S_zz = result.seebeck[3, 3, iT, iμ] * 1e6

                κ_xx = result.kappa[1, 1, iT, iμ]
                κ_yy = result.kappa[2, 2, iT, iμ]
                κ_zz = result.kappa[3, 3, iT, iμ]

                println(io, "$T,$μ,$σ_xx,$σ_yy,$σ_zz,$S_xx,$S_yy,$S_zz,$κ_xx,$κ_yy,$κ_zz")
            end
        end
    end
end

# ============================================================================
# Legacy API (for backward compatibility)
# ============================================================================

#=
    save_calculation(filename, data, equivalences, coeffs; metadata=nothing)

Legacy API - saves interpolation data.
Prefer using `save_interpolation` with InterpolationResult or FourierInterpolator.
=#
function save_calculation(
    filename::String,
    lattvec,
    equivalences,
    coeffs;
    atoms = nothing,
    metadata = nothing,
)
    meta = isnothing(metadata) ? Dict{String,Any}() : metadata
    equivs = [Matrix{Int}(eq) for eq in equivalences]
    result = InterpolationResult(coeffs, equivs, Matrix{Float64}(lattvec), atoms, meta)
    save_interpolation(filename, result)
end

#=
    load_calculation(filename) -> (lattvec, equivalences, coeffs, atoms, metadata)

Legacy API - loads interpolation data as tuple.
Prefer using `load_interpolation` which returns InterpolationResult.
=#
function load_calculation(filename::String)
    result = load_interpolation(filename)
    return (
        result.lattvec,
        result.equivalences,
        result.coeffs,
        result.atoms,
        result.metadata,
    )
end
