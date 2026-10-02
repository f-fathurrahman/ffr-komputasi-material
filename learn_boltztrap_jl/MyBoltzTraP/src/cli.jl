# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using Comonicon: @cast, @main
using JLD2: load
using Logging

#=
    _setup_logging(debug::Bool)

Configure logging based on --debug flag.
If debug=true, enable Debug level logging for BoltzTraP module.
=#
function _setup_logging(debug::Bool)
    if debug
        # Enable debug logging
        ENV["JULIA_DEBUG"] = "BoltzTraP"
        global_logger(ConsoleLogger(stderr, Logging.Debug))
    end
end

"""
    interpolate(directory; format="auto", output="", kpoints=0, multiplier=0, emin=-Inf, emax=Inf, absolute=false, verbose=false, debug=false)

Interpolate band structure from DFT data.

Read DFT output and compute Fourier interpolation coefficients
for band energies. Support VASP, Quantum ESPRESSO, and other formats.

# Arguments

  - `directory`: Directory containing DFT output files

# Options

  - `-f, --format <fmt>`: DFT format (auto, vasp, qe). Default: auto-detect
  - `-o, --output <file>`: Output file path (default: based on input directory)
  - `-k, --kpoints <n>`: Target number of k-points/equivalences (default: 5000)
  - `-m, --multiplier <n>`: Multiplier for k-points (alternative to --kpoints)
  - `--emin <e>`: Minimum energy relative to Fermi in Ha (default: -Inf)
  - `--emax <e>`: Maximum energy relative to Fermi in Ha (default: +Inf)

# Flags

  - `--absolute`: Interpret emin/emax as absolute energies
  - `-v, --verbose`: Print progress information
  - `--debug`: Enable debug logging (detailed internal info)

# Examples

```bash
boltztrap interpolate ./Si.vasp
boltztrap interpolate ./Si.vasp -f vasp -o si_interp.jld2 -k 10000
boltztrap interpolate ./Si.qe -f qe --emin -0.5 --emax 0.5 -v
boltztrap interpolate ./calculation --debug  # with debug output
```
"""
@cast function interpolate(
    directory::String;
    format::String = "auto",
    output::String = "",
    kpoints::Int = 0,
    multiplier::Int = 0,
    emin::Float64 = -Inf,
    emax::Float64 = Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    debug::Bool = false,
)
    _setup_logging(debug)
    @debug "CLI interpolate called" directory format output kpoints multiplier emin emax absolute verbose debug

    # Determine output filename
    if isempty(output)
        base = basename(rstrip(directory, '/'))
        output = base * "_interp.jld2"
    end

    # Determine kpoints setting
    kpts = kpoints > 0 ? kpoints : nothing
    mult = multiplier > 0 ? multiplier : nothing

    if isnothing(kpts) && isnothing(mult)
        kpts = 5000  # Default
    end

    verbose && println("BoltzTraP.jl interpolate")
    verbose && println("  Input: $directory")
    verbose && println("  Format: $format")
    verbose && println("  Output: $output")

    # Load DFT data based on format
    format_lower = lowercase(format)
    data = if format_lower == "auto"
        verbose && println("  Auto-detecting format...")
        detected = detected_format(directory)
        if isnothing(detected)
            error("Could not auto-detect DFT format in $directory")
        end
        verbose && println("  Detected: $detected")
        load_dft(directory)
    elseif format_lower == "vasp"
        load_vasp(directory)
    elseif format_lower == "qe"
        load_qe(directory)
    else
        error("Unknown format: $format. Supported: auto, vasp, qe")
    end

    @debug "DFT data loaded" format nkpts=size(data.kpoints, 2) nbands=size(data.ebands, 1)

    # Run interpolation
    result = run_interpolate(
        data;
        source = directory,
        output = output,
        kpoints = kpts,
        multiplier = mult,
        emin = emin,
        emax = emax,
        absolute = absolute,
        verbose = verbose,
    )

    println("Saved interpolation to: $output")
    println("  Equivalences: $(length(result.equivalences))")
    println("  Bands: $(size(result.coeffs, 1))")

    return nothing
end

"""
    integrate(input; temperature="300", output="", bins=0, scissor=NaN, verbose=false, debug=false)

Compute transport coefficients from interpolated band structure.

Read interpolation result and compute electrical conductivity,
Seebeck coefficient, and thermal conductivity.

# Arguments

  - `input`: Path to interpolation result file (.jld2 or .bt2)

# Options

  - `-t, --temperature <T>`: Temperature(s) in K. Formats:

      + Single value: `"300"` → [300.0]
      + Range: `"100:500:50"` (start:stop:step) → [100.0, 150.0, 200.0, ..., 500.0]
      + List: `"100,200,300"` → [100.0, 200.0, 300.0]

  - `-o, --output <file>`: Output file path (default: based on input)

  - `-b, --bins <n>`: Number of DOS histogram bins (default: auto)

  - `-s, --scissor <eV>`: Target band gap in eV for scissor correction.
    Shifts conduction bands to achieve the specified gap. Only for semiconductors.

# Flags

  - `-v, --verbose`: Print progress information
  - `--debug`: Enable debug logging (detailed internal info)

# Examples

```bash
boltztrap integrate si_interp.jld2 -t 300
boltztrap integrate si_interp.jld2 -t "100:500:50" -o results.jld2
boltztrap integrate si_interp.jld2 -t 300 --scissor 1.1
boltztrap integrate si_interp.jld2 -t "100,200,300,400,500" -v --debug
```
"""
@cast function integrate(
    input::String;
    temperature::String = "300",
    output::String = "",
    bins::Int = 0,
    scissor::Float64 = NaN,
    verbose::Bool = false,
    debug::Bool = false,
)
    _setup_logging(debug)
    @debug "CLI integrate called" input temperature output bins scissor verbose debug

    # Parse temperature string
    temperatures = _parse_temperatures(temperature)

    # Determine output filename
    if isempty(output)
        base = splitext(basename(input))[1]
        output = base * "_transport.jld2"
    end

    # Convert NaN to nothing for scissor
    scissor_arg = isnan(scissor) ? nothing : scissor

    verbose && println("BoltzTraP.jl integrate")
    verbose && println("  Input: $input")
    verbose && println("  Output: $output")
    verbose && println("  Temperatures: $temperatures K")
    if !isnothing(scissor_arg)
        verbose && println("  Scissor: $scissor_arg eV")
    end

    # Run integration
    result = run_integrate(
        input;
        temperatures = temperatures,
        output = output,
        bins = bins,
        scissor = scissor_arg,
        verbose = verbose,
    )

    println("Saved transport results to: $output")
    println("  Temperatures: $(length(result.temperatures))")
    println("  μ points: $(length(result.mu_values))")

    return nothing
end

#=
    _parse_temperatures(s::String) -> Vector{Float64}

Parse temperature string into vector of temperatures.

Supported formats:
- Single value: "300" -> [300.0]
- Range: "100:100:500" -> [100.0, 200.0, ..., 500.0] (Julia range notation: start:step:stop)
- List: "100,200,300" -> [100.0, 200.0, 300.0]
=#
function _parse_temperatures(s::String)::Vector{Float64}
    s = strip(s)

    if occursin(':', s)
        # Range format: start:step:stop (Julia notation)
        parts = split(s, ':')
        if length(parts) == 3
            start = parse(Float64, parts[1])
            step = parse(Float64, parts[2])
            stop = parse(Float64, parts[3])
            return collect(start:step:stop)
        elseif length(parts) == 2
            start = parse(Float64, parts[1])
            stop = parse(Float64, parts[2])
            return collect(start:100.0:stop)  # Default step
        else
            error("Invalid temperature range format: $s")
        end
    elseif occursin(',', s)
        # List format: T1,T2,T3
        parts = split(s, ',')
        return [parse(Float64, strip(p)) for p in parts]
    else
        # Single value
        return [parse(Float64, s)]
    end
end

"""
    describe(file)

Display information about a BoltzTraP.jl result file.

Show basic information and metadata from interpolation (.jld2/.bt2) or
transport result files (.jld2).

# Arguments

  - `file`: Path to result file (.jld2 or .bt2)

# Examples

```bash
boltztrap describe si_interp.jld2
boltztrap describe si_interp.bt2
boltztrap describe si_transport.jld2
```
"""
@cast function describe(file::String)
    # Delegate to internal function (API describe works on result objects)
    _describe_file(file)
    return nothing
end

#=
    _parse_component(s::String) -> Tuple{Int, Int}

Parse tensor component string to indices.
"xx" -> (1,1), "yy" -> (2,2), "zz" -> (3,3)
"xy" -> (1,2), "xz" -> (1,3), "yz" -> (2,3)
=#
function _parse_component(s::String)::Tuple{Int,Int}
    s = lowercase(strip(s))
    component_map = Dict(
        "xx" => (1, 1),
        "yy" => (2, 2),
        "zz" => (3, 3),
        "xy" => (1, 2),
        "yx" => (2, 1),
        "xz" => (1, 3),
        "zx" => (3, 1),
        "yz" => (2, 3),
        "zy" => (3, 2),
    )
    if haskey(component_map, s)
        return component_map[s]
    else
        error("Invalid component: $s. Valid: xx, yy, zz, xy, xz, yz")
    end
end

#=
    _get_quantity_data(data::Dict, quantity::String) -> Array

Extract quantity array from transport result data.
=#
function _get_quantity_data(data::Dict, quantity::String)
    quantity = lowercase(strip(quantity))
    quantity_map = Dict(
        "sigma" => "sigma",
        "conductivity" => "sigma",
        "seebeck" => "seebeck",
        "s" => "seebeck",
        "kappa" => "kappa",
        "thermal" => "kappa",
        "dos" => "dos",
        "n" => "N",
        "carrier" => "N",
    )

    key = get(quantity_map, quantity, quantity)
    if !haskey(data, key)
        available = filter(k -> k in keys(data), ["sigma", "seebeck", "kappa", "dos", "N"])
        error("Quantity '$quantity' not found. Available: $(join(available, ", "))")
    end
    return data[key]
end

#=
    _get_quantity_label(quantity::String) -> String

Get display label for quantity with units.
=#
function _get_quantity_label(quantity::String)
    quantity = lowercase(strip(quantity))
    labels = Dict(
        "sigma" => "σ/τ (S/m)",
        "seebeck" => "S (μV/K)",
        "kappa" => "κ/τ (W/m/K)",
        "dos" => "DOS (states/Ha/cell)",
        "n" => "n (carriers/cell)",
    )
    return get(labels, quantity, quantity)
end

"""
    plot(file; quantity="sigma", component="xx", abscissa="mu", temperature=300.0, mu=NaN, output="")

Plot transport coefficients from integration result.

Plot transport properties (σ, S, κ, etc.) as a function of
chemical potential or temperature.

# Arguments

  - `file`: Path to transport result file (.jld2)

# Options

  - `-q, --quantity <name>`: Quantity to plot (sigma, seebeck, kappa, dos, n). Default: sigma
  - `-c, --component <ij>`: Tensor component (xx, yy, zz, xy, xz, yz). Default: xx
  - `-a, --abscissa <var>`: X-axis variable (mu or T). Default: mu
  - `-t, --temperature <T>`: Temperature in K for mu plot. Default: 300
  - `-m, --mu <μ>`: Chemical potential in Ha for T plot. Default: Fermi energy
  - `-o, --output <file>`: Output file (PNG/PDF). Default: display

# Examples

```bash
boltztrap plot transport.jld2 -q sigma -c xx
boltztrap plot transport.jld2 -q seebeck -a T -m 0.0
boltztrap plot transport.jld2 -q sigma -t 300 -o sigma_300K.png
```
"""
@cast function plot(
    file::String;
    quantity::String = "sigma",
    component::String = "xx",
    abscissa::String = "mu",
    temperature::Float64 = 300.0,
    mu::Float64 = NaN,
    output::String = "",
)
    if !isfile(file)
        error("File not found: $file")
    end

    if !endswith(file, ".jld2")
        error("Expected .jld2 file, got: $file")
    end

    # Load data
    data = load(file)

    # Verify it's a transport result
    if !haskey(data, "temperatures") || !haskey(data, "sigma")
        error("Not a transport result file. Use 'describe' to check file type.")
    end

    temperatures = data["temperatures"]
    mu_values = data["mu_values"]

    # Get quantity data
    qdata = _get_quantity_data(data, quantity)
    qlabel = _get_quantity_label(quantity)

    # Parse component
    i, j = _parse_component(component)

    # Get Fermi energy from metadata if available
    # Prefer fermi_dft_Ha (original DFT Fermi) over fermi_Ha (refined Fermi)
    fermi = 0.0
    if haskey(data, "metadata")
        md = data["metadata"]
        if haskey(md, "fermi_dft_Ha")
            fermi = md["fermi_dft_Ha"]
        elseif haskey(md, "fermi_Ha")
            fermi = md["fermi_Ha"]
        elseif haskey(md, "fermi_energy")
            fermi = md["fermi_energy"]
        end
    end

    abscissa = lowercase(strip(abscissa))

    if abscissa == "mu"
        # Plot vs chemical potential at fixed temperature
        # Find closest temperature
        t_idx = argmin(abs.(temperatures .- temperature))
        actual_T = temperatures[t_idx]

        # Extract data: qdata has shape (3, 3, nT, nmu) for tensors
        # or (nT, nmu) for scalars like DOS
        if ndims(qdata) == 4
            y_data = qdata[i, j, t_idx, :]
        elseif ndims(qdata) == 2
            y_data = qdata[t_idx, :]
        else
            error("Unexpected data shape: $(size(qdata))")
        end

        # Convert Seebeck to μV/K
        if lowercase(quantity) in ["seebeck", "s"]
            y_data = y_data .* 1e6  # V/K → μV/K
        end

        # mu_values is already in eV, fermi is in Ha
        # Convert fermi to eV and compute relative μ (HA_TO_EV from units.jl)
        fermi_eV = fermi * HA_TO_EV
        mu_eV = mu_values .- fermi_eV  # Both in eV now

        p = Plots.plot(
            mu_eV,
            y_data;
            xlabel = "μ - εF (eV)",
            ylabel = qlabel,
            title = "$(uppercase(quantity)) at T = $(round(actual_T; digits=1)) K",
            legend = false,
            linewidth = 2,
        )

        # Add vertical line at Fermi level
        Plots.vline!(p, [0.0]; linestyle = :dash, color = :gray, alpha = 0.5)

    elseif abscissa == "t"
        # Plot vs temperature at fixed chemical potential
        # mu_values is in eV, fermi is in Ha (HA_TO_EV from units.jl)
        fermi_eV = fermi * HA_TO_EV

        # Use provided mu (in eV) or Fermi energy (converted to eV)
        target_mu = isnan(mu) ? fermi_eV : mu

        # Find closest mu value (both in eV)
        mu_idx = argmin(abs.(mu_values .- target_mu))
        actual_mu = mu_values[mu_idx]

        # Extract data: qdata has shape (3, 3, nT, nmu) for tensors
        if ndims(qdata) == 4
            y_data = qdata[i, j, :, mu_idx]
        elseif ndims(qdata) == 2
            y_data = qdata[:, mu_idx]
        else
            error("Unexpected data shape: $(size(qdata))")
        end

        # Convert Seebeck to μV/K
        if lowercase(quantity) in ["seebeck", "s"]
            y_data = y_data .* 1e6  # V/K → μV/K
        end

        mu_rel_eV = actual_mu - fermi_eV  # Both in eV

        p = Plots.plot(
            temperatures,
            y_data;
            xlabel = "Temperature (K)",
            ylabel = qlabel,
            title = "$(uppercase(quantity)) at μ - εF = $(round(mu_rel_eV; digits=3)) eV",
            legend = false,
            linewidth = 2,
        )
    else
        error("Invalid abscissa: $abscissa. Use 'mu' or 'T'")
    end

    # Save or display
    if !isempty(output)
        Plots.savefig(p, output)
        println("Saved plot to: $output")
    else
        display(p)
    end

    return nothing
end

#=
    _parse_kpath_string(kpath_str) -> NamedTuple

Parse k-path string format: "G:0,0,0;X:0.5,0,0.5;L:0.5,0.5,0.5|G-X-L-G"
- Points section (before |): label:x,y,z;...
- Paths section (after |): label-label-...|label-label-...
=#
function _parse_kpath_string(kpath_str::String)
    parts = split(kpath_str, "|")
    if length(parts) < 2
        error(
            "Invalid kpath format. Expected: 'label:x,y,z;...|path1-path2-...' " *
            "Example: 'G:0,0,0;X:0.5,0,0.5;L:0.5,0.5,0.5|G-X-L-G'",
        )
    end

    # Parse points
    points_str = parts[1]
    points = Dict{String,Vector{Float64}}()
    for pt in split(points_str, ";")
        pt = strip(pt)
        isempty(pt) && continue
        label_coords = split(pt, ":")
        if length(label_coords) != 2
            error("Invalid point format: '$pt'. Expected 'label:x,y,z'")
        end
        label = strip(label_coords[1])
        coords = parse.(Float64, split(label_coords[2], ","))
        if length(coords) != 3
            error("Invalid coordinates for '$label': expected 3 values")
        end
        points[label] = coords
    end

    # Parse paths (remaining parts)
    path_segments = Vector{Vector{String}}()
    for i = 2:length(parts)
        path_str = strip(parts[i])
        isempty(path_str) && continue
        labels = String.(strip.(split(path_str, "-")))
        # Validate labels exist
        for lbl in labels
            if !haskey(points, lbl)
                error("Unknown point label '$lbl' in path. Define it in points section.")
            end
        end
        push!(path_segments, labels)
    end

    if isempty(path_segments)
        error("No path segments found in kpath string")
    end

    return (points = points, paths = path_segments)
end

"""
    plotbands(file; npoints=100, emin=-5.0, emax=5.0, no_fermi=false, kpath="", output="")

Plot band structure along high-symmetry k-path.

Generate band structure plot from interpolation result using
Brillouin.jl to determine the k-path for the crystal's space group.
For .bt2 files without spacegroup metadata, use --kpath to specify manually.

# Arguments

  - `file`: Path to interpolation result file (.jld2 or .bt2)

# Options

  - `-n, --npoints <n>`: Number of k-points per path segment. Default: 100
  - `--emin <e>`: Minimum energy relative to Fermi [eV]. Default: -5.0
  - `--emax <e>`: Maximum energy relative to Fermi [eV]. Default: 5.0
  - `--kpath <spec>`: Manual k-path specification (required for .bt2 without spacegroup)
  - `-o, --output <file>`: Output file (PNG/PDF). Default: display

# Flags

  - `--no-fermi`: Hide Fermi level line

# K-path format

`--kpath "label:x,y,z;...|path-path-..."`

Example for FCC (Si): `--kpath "G:0,0,0;X:0.5,0,0.5;L:0.5,0.5,0.5;W:0.5,0.25,0.75|G-X-W-L-G"`

# Examples

```bash
boltztrap plotbands si_interp.jld2
boltztrap plotbands si_interp.jld2 --emin -10 --emax 10 -o bands.png
boltztrap plotbands si_interp.jld2 -n 200 --no-fermi
boltztrap plotbands si_interp.bt2 --kpath "G:0,0,0;X:0.5,0,0.5|G-X-G"
```
"""
@cast function plotbands(
    file::String;
    npoints::Int = 100,
    emin::Float64 = -5.0,
    emax::Float64 = 5.0,
    no_fermi::Bool = false,
    kpath::String = "",
    output::String = "",
)
    # Load Plots and Brillouin for plotting extension
    try
        @eval using Plots
        @eval using Brillouin
    catch
        error(
            "plotbands requires Plots.jl and Brillouin.jl. Install with: Pkg.add([\"Plots\", \"Brillouin\"])",
        )
    end

    if !isfile(file)
        error("File not found: $file")
    end

    ext = splitext(file)[2]
    if ext ∉ (".jld2", ".bt2")
        error("Unsupported file format: $ext. Expected .jld2 or .bt2")
    end

    # Load interpolation result
    result = load_interpolation(file)

    # Determine k-path source
    kpath_arg = nothing
    if !isempty(kpath)
        # Manual k-path specified
        kpath_arg = _parse_kpath_string(kpath)
        println("Plotting band structure (manual k-path)...")
        println("  Points: $(join(keys(kpath_arg.points), ", "))")
        println("  Paths: $(join([join(p, "-") for p in kpath_arg.paths], " | "))")
    else
        # Auto k-path from spacegroup
        if !haskey(result.metadata, "spacegroup_number")
            error(
                "Space group not found in file. " *
                "Use --kpath to specify k-path manually, or use a .jld2 file with spacegroup metadata.\n" *
                "Example: --kpath \"G:0,0,0;X:0.5,0,0.5;L:0.5,0.5,0.5|G-X-L-G\"",
            )
        end
        println("Plotting band structure...")
        println(
            "  Space group: $(result.metadata["spacegroup_number"]) ($(get(result.metadata, "spacegroup_symbol", "")))",
        )
    end

    println("  K-points per segment: $npoints")
    println("  Energy range: [$emin, $emax] eV")

    # Generate plot
    out = isempty(output) ? nothing : output
    p = plot_bands(
        result;
        npoints = npoints,
        emin = emin,
        emax = emax,
        fermi_line = !no_fermi,
        kpath = kpath_arg,
        output = out,
    )

    # Display if no output file
    if isempty(output)
        display(p)
    end

    return nothing
end

"""
BoltzTraP.jl - Band structure interpolation and transport coefficients.

Julia implementation of BoltzTraP2.

# Commands

  - `interpolate`: Interpolate band structure from DFT data
  - `integrate`: Compute transport coefficients
  - `describe`: Display file information
  - `plot`: Plot transport coefficients
  - `plotbands`: Plot band structure along k-path

# Examples

```bash
boltztrap interpolate ./Si.vasp -o si.jld2
boltztrap integrate si.jld2 -t 300
boltztrap describe si_interp.jld2
boltztrap plot si_transport.jld2 -q sigma -c xx
boltztrap plotbands si_interp.jld2 -o bands.png
boltztrap --help
```
"""
@main
