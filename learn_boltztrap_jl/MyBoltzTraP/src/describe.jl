# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using JLD2

"""
    describe(result::InterpolationResult; io::IO=stdout)

Display summary of interpolation result.

# Arguments

  - `result`: [`InterpolationResult`](@ref) from [`run_interpolate`](@ref) or [`load_interpolation`](@ref)
  - `io=stdout`: Output stream

# Example

```julia
interp = run_interpolate(\"./Si.vasp\")
describe(interp)
```
"""
function describe(result::InterpolationResult; io::IO = stdout)
    println(io, "=" ^ 50)
    println(io, "Interpolation Result")
    println(io, "=" ^ 50)
    println(io)

    # Basic info
    println(io, "Data:")
    println(io, "  Bands: $(size(result.coeffs, 1))")
    println(io, "  Equivalences: $(length(result.equivalences))")
    println(io, "  Lattice vectors: $(size(result.lattvec))")
    println(io)

    # Metadata
    if !isempty(result.metadata)
        println(io, "Metadata:")
        for (k, v) in sort(collect(result.metadata), by = first)
            v_str = _format_value(v)
            println(io, "  $k: $v_str")
        end
    end

    return nothing
end

"""
    describe(result::TransportResult; io::IO=stdout)

Display summary of transport result.

# Arguments

  - `result`: [`TransportResult`](@ref) from [`run_integrate`](@ref) or [`load_integrate`](@ref)
  - `io=stdout`: Output stream

# Example

```julia
transport = run_integrate(interp; temperatures = [300.0])
describe(transport)
```
"""
function describe(result::TransportResult; io::IO = stdout)
    println(io, "=" ^ 50)
    println(io, "Transport Result")
    println(io, "=" ^ 50)
    println(io)

    # Basic info
    println(io, "Data:")
    println(io, "  Temperatures: $(length(result.temperatures)) points")
    if length(result.temperatures) <= 5
        println(io, "    $(result.temperatures) K")
    else
        println(io, "    $(result.temperatures[1]) - $(result.temperatures[end]) K")
    end
    println(io, "  Chemical potentials: $(length(result.mu_values)) points")
    println(io, "    $(minimum(result.mu_values)) to $(maximum(result.mu_values)) eV")
    println(io, "  σ tensor shape: $(size(result.sigma))")
    println(io)

    # Available quantities
    quantities = String[]
    !isempty(result.seebeck) && push!(quantities, "seebeck")
    !isempty(result.kappa) && push!(quantities, "kappa")
    !isnothing(result.dos) && push!(quantities, "dos")

    if !isempty(quantities)
        println(io, "Available quantities:")
        for q in quantities
            println(io, "  - $q")
        end
        println(io)
    end

    # Metadata
    if !isempty(result.metadata)
        println(io, "Metadata:")
        for (k, v) in sort(collect(result.metadata), by = first)
            v_str = _format_value(v)
            println(io, "  $k: $v_str")
        end
    end

    return nothing
end

#=
    _describe_file(file::String; io::IO=stdout)

Internal function to load and describe a result file (.jld2 or .bt2).
Called by CLI describe command.
=#
function _describe_file(file::String; io::IO = stdout)
    if !isfile(file)
        error("File not found: $file")
    end

    ext = splitext(file)[2]

    if ext == ".bt2"
        # Python BoltzTraP2 format - always interpolation result
        result = load_interpolation(file)
        println(io, "File: $file")
        println(io, "Format: Python BoltzTraP2 (.bt2)")
        println(io)
        describe(result; io = io)
    elseif ext == ".jld2"
        # Julia format - detect type by trying to load
        # First try interpolation
        try
            result = load_interpolation(file)
            println(io, "File: $file")
            println(io, "Format: JLD2 (.jld2)")
            println(io)
            describe(result; io = io)
            return nothing
        catch
        end

        # Then try transport
        try
            result = load_integrate(file)
            println(io, "File: $file")
            println(io, "Format: JLD2 (.jld2)")
            println(io)
            describe(result; io = io)
            return nothing
        catch
        end

        # Unknown format - show raw keys
        data = JLD2.load(file)
        println(io, "File: $file")
        println(io, "Type: Unknown")
        println(io, "Keys: $(collect(keys(data)))")
    else
        error("Unsupported file format: $ext. Expected .jld2 or .bt2")
    end

    return nothing
end #= Helper function to format values for display =#

function _format_value(v)
    if v isa AbstractFloat
        string(round(v; sigdigits = 6))
    elseif v isa AbstractRange
        "$(first(v)):$(last(v))"
    else
        string(v)
    end
end
