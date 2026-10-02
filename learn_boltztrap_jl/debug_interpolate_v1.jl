using MyBoltzTraP

# Some internal functions not exported
using MyBoltzTraP: get_spacegroup_info, get_equivalences

function my_run_interpolate(
    data::DFTData{1};
    source::String = "unknown",
    output::Union{String,Nothing} = nothing,
    kpoints::Union{Int,Nothing} = nothing,
    multiplier::Union{Int,Nothing} = nothing,
    emin::Float64 = -Inf,
    emax::Float64 = +Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    symprec::Real = 1e-5,
    dosweight::Union{Float64,Nothing} = nothing,
)
    #ffr: really need to convert to NamedTuple ?
    # Convert DFTData to NamedTuple for existing implementation
    data_nt = (
        lattice = data.lattice,
        positions = data.positions,
        species = data.species,
        kpoints = data.kpoints,
        weights = data.weights,
        ebands = data.ebands,
        occupations = data.occupations,
        fermi = data.fermi,
        nelect = data.nelect,
    )
    @info "Pass here 102"
    return my_run_interpolate(
        data_nt;
        source,
        output,
        kpoints,
        multiplier,
        emin,
        emax,
        absolute,
        verbose,
        symprec,
        dosweight,
    )
end

# This is the actual driver for run_interpolate.
# Main input is a NamedTuple.
function my_run_interpolate(
    data::NamedTuple;
    source::String = "unknown",
    output::Union{String,Nothing} = nothing,
    kpoints::Union{Int,Nothing} = nothing,
    multiplier::Union{Int,Nothing} = nothing,
    emin::Float64 = -Inf,
    emax::Float64 = +Inf,
    absolute::Bool = false,
    verbose::Bool = false,
    symprec::Real = 1e-5,
    dosweight::Union{Float64,Nothing} = nothing,
)

    println("\n<div> ENTER my_run_interpolate\n")

    # Log received arguments for debugging
    @info "run_interpolate called" source output kpoints multiplier emin emax absolute verbose symprec

    # Validate arguments
    if isnothing(kpoints) && isnothing(multiplier)
        kpoints = 5000  # Default
    elseif !isnothing(kpoints) && !isnothing(multiplier)
        error("Cannot specify both kpoints and multiplier")
    end

    @info "run_interpolate after defaults" kpoints multiplier

    verbose && println("Processing data from $source...")

    # 1. Prepare structure data
    lattvec = data.lattice
    positions = data.positions  # 3×natoms, need to transpose
    species = data.species
    types = [findfirst(==(s), unique(species)) for s in species]
    magmom = get(data, :magmom, nothing)  # For collinear magnetic calculations

    # 2. Determine target k-points
    if !isnothing(multiplier)
        nkpt_target = size(data.kpoints, 2) * multiplier
        verbose && println("Using multiplier=$multiplier → nkpt_target=$nkpt_target")
    else
        nkpt_target = kpoints
    end

    # 3. Get space group info
    sginfo = get_spacegroup_info(lattvec, positions', types; symprec)
    verbose && println(
        "  Space group: $(sginfo.spacegroup_number) ($(sginfo.international_symbol))",
    )

    # 4. Compute equivalences
    verbose && println("Computing equivalences for ~$nkpt_target k-points...")
    equivalences, radius, nrot =
        get_equivalences(lattvec, positions', types, magmom, nkpt_target; symprec)
    verbose && println("  Rotations: $nrot")
    verbose && println("  Radius: $(round(radius, digits=2))")
    verbose && println("  Equivalences: $(length(equivalences))")

    # 4. Prepare band data (already in atomic units)
    kpts = data.kpoints'  # nk×3
    ebands_raw = data.ebands[:, :, 1]  # nbands×nk (spin 1), in Ha
    fermi = data.fermi  # in Ha
    nbands_total = size(ebands_raw, 1)

    # Determine dosweight from spin polarization (or use explicit override)
    # Use magmom to detect collinear (ebands may be concatenated for collinear)
    # Override with explicit dosweight for SOC (dosweight=1.0, magmom=nothing)
    if isnothing(dosweight)
        dosweight = isnothing(magmom) ? 2.0 : 1.0
    end

    # 5. Filter bands by energy
    # emin/emax are in Ha (relative to Fermi unless absolute=true)
    if absolute
        emin_abs = emin
        emax_abs = emax
    else
        emin_abs = fermi + emin
        emax_abs = fermi + emax
    end

    # Find bands within energy range
    band_min = minimum(ebands_raw, dims = 2)[:]
    band_max = maximum(ebands_raw, dims = 2)[:]
    band_mask = (band_max .>= emin_abs) .& (band_min .<= emax_abs)
    selected_bands = findall(band_mask)

    if isempty(selected_bands)
        error("No bands found in energy range [$emin_abs, $emax_abs] Ha")
    end

    ebands = ebands_raw[selected_bands, :]
    verbose && println("  Bands: $(length(selected_bands))/$nbands_total selected")
    verbose && println(
        "  Energy range: [$(round(emin_abs, digits=6)), $(round(emax_abs, digits=6))] Ha",
    )

    # 6. Convert equivalences to matrix format for FourierInterpolator
    equiv_matrices = [hcat(eq...)' for eq in equivalences]

    # 7. Create interpolator
    verbose && println("Fitting Fourier coefficients...")
    interp = FourierInterpolator(kpts, ebands, equiv_matrices, lattvec)

    # 8. Create result with metadata
    atoms = Dict{String,Any}(
        "species" => species,
        "positions" => collect(positions'),
        "lattice" => collect(lattvec),
    )
    # Determine spintype from magmom
    spintype = isnothing(magmom) ? "Unpolarized" : "Collinear"

    metadata = Dict{String,Any}(
        "fermi" => fermi,
        "nelect" => data.nelect,
        "dosweight" => dosweight,
        "spintype" => spintype,
        "selected_bands" => selected_bands,
        "nkpt_original" => size(data.kpoints, 2),
        "nkpt_target" => nkpt_target,
        "radius" => radius,
        "nrotations" => nrot,
        "emin" => emin,
        "emax" => emax,
        "absolute" => absolute,
        "source" => source,
        "spacegroup_number" => sginfo.spacegroup_number,
        "spacegroup_symbol" => sginfo.international_symbol,
    )

    result = InterpolationResult(interp; atoms = atoms, metadata = metadata)

    # 9. Save if output specified
    if !isnothing(output)
        verbose && println("Saving to $output...")
        save_interpolation(output, result)
    end

    verbose && println("Done.")

    println("\n</div> EXIT my_run_interpolate\n")

    return result
end


function debug_main()
    # Path to Si VASP data (directory containing vasprun.xml and POSCAR)
    datadir = "./data_Si_vasp"
    data = load_vasp(datadir)
    println("Type of data = ", typeof(data))
    interp = my_run_interpolate(data; kpoints = 5000, verbose = true)
    return
end

