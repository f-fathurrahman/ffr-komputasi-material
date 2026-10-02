# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Logging

#=
    LoaderError <: Exception

Exception thrown when a loader fails to read data from a directory.
=#
struct LoaderError <: Exception
    msg::String
    loader::String
end

function Base.showerror(io::IO, e::LoaderError)
    print(io, "LoaderError in $(e.loader): $(e.msg)")
end

#=
    Loader registry

Each loader is a tuple of (name, detector, loader_func):
- name: Format name (e.g., "VASP", "QE")
- detector: Function that returns `true` if the format is detected
- loader_func: Function that loads the data

Loaders are tried in reverse order (last registered = highest priority).
=#
const LOADERS = Tuple{String,Function,Function}[]

#=
    register_loader(name::String, detector::Function, loader::Function)

Register a new DFT data loader.

# Arguments
- `name`: Format name (e.g., "VASP", "QE")
- `detector`: Function `(directory) -> Bool` that checks if format is present
- `loader`: Function `(directory) -> DFTData` that loads the data
=#
function register_loader(name::String, detector::Function, loader::Function)
    push!(LOADERS, (name, detector, loader))
    @debug "Registered loader" name
end

# Detector functions
function _detect_vasp(directory::String)
    isfile(joinpath(directory, "vasprun.xml"))
end

function _detect_qe(directory::String)
    # QE uses data-file-schema.xml or data-file.xml in *.save directory
    for item in readdir(directory; join = true)
        if isdir(item) && endswith(item, ".save")
            if isfile(joinpath(item, "data-file-schema.xml")) ||
               isfile(joinpath(item, "data-file.xml"))
                return true
            end
        end
    end
    # Also check for direct XML files
    for item in readdir(directory; join = true)
        if isfile(item) && endswith(item, ".xml")
            # Quick check for QE signature in first 10 lines
            try
                found = false
                open(item) do f
                    for _ = 1:10
                        eof(f) && break
                        line = readline(f)
                        if occursin("qes:", line) || occursin("espresso", lowercase(line))
                            found = true
                            break
                        end
                    end
                end
                if found
                    return true
                end
            catch
            end
        end
    end
    return false
end

# Register built-in loaders (order matters: last = highest priority)
# GENE is generic format, so it has lowest priority
function _init_loaders()
    if isempty(LOADERS)
        register_loader("GENE", _detect_gene, load_gene)
        register_loader("ABINIT", _detect_abinit, load_abinit)
        register_loader("Wien2k", _detect_wien2k, load_wien2k)
        register_loader("QE", _detect_qe, load_qe)
        register_loader("VASP", _detect_vasp, load_vasp)
    end
end

"""
    load_dft(directory::String) -> [`DFTData`](@ref)

Auto-detect DFT format and load band structure data.

Try registered loaders in reverse order (VASP first, then QE, etc.).
Throw an error if no compatible format is found.

# Arguments

  - `directory`: Path to directory containing DFT output files

# Returns

  - [`DFTData`](@ref) with band structure data (same format as [`load_vasp`](@ref), [`load_qe`](@ref)).

# Example

```julia
# Auto-detect format
data = load_dft(\"./calculation\")

# Explicit format (still works)
data = load_vasp(\"./Si.vasp\")
data = load_qe(\"./Si.qe\")
```
"""
function load_dft(directory::String)
    _init_loaders()

    if !isdir(directory)
        error("Directory not found: $directory")
    end

    @debug "load_dft: trying loaders" directory n_loaders=length(LOADERS)

    errors = String[]

    # Try loaders in reverse order (last registered = highest priority)
    for (name, detector, loader) in reverse(LOADERS)
        @debug "Trying loader" name
        try
            if detector(directory)
                @debug "Format detected, loading" name
                data = loader(directory)
                @info "Successfully loaded $name data from $directory"
                return data
            else
                @debug "Format not detected" name
            end
        catch e
            msg = e isa LoaderError ? e.msg : string(e)
            push!(errors, "$name: $msg")
            @debug "Loader failed" name error=msg
        end
    end

    # No loader succeeded
    error_details = isempty(errors) ? "" : "\nLoader errors:\n  " * join(errors, "\n  ")
    error("No compatible DFT format found in $directory$error_details")
end

"""
    detected_format(directory::String) -> Union{String, Nothing}

Check which DFT format is detected in the directory without loading.

# Returns

Format name (e.g., "VASP", "QE") or `nothing` if no format detected.
"""
function detected_format(directory::String)
    _init_loaders()

    if !isdir(directory)
        return nothing
    end

    for (name, detector, _) in reverse(LOADERS)
        try
            if detector(directory)
                return name
            end
        catch
        end
    end
    return nothing
end
