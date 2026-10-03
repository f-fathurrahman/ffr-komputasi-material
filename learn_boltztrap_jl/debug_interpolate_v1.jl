using StaticArrays
import Dates
using LinearAlgebra: norm, det, pinv
using Serialization: serialize, deserialize
import Spglib

# These are the functions/types that are still used.
using MyBoltzTraP: SpinType, AbstractInterpolator, DFTData
using MyBoltzTraP: Unpolarized, Collinear, NonCollinear
using MyBoltzTraP: load_vasp

function my_get_spintype(magmom)
    if isnothing(magmom)
        return Unpolarized()
    elseif eltype(magmom) <: Real
        return Collinear()
    else
        return NonCollinear()
    end
end
# Backward compatibility alias
const my_get_magmom_type = my_get_spintype

function my_determine_compatibility(
    rotation::AbstractMatrix,
    perm::AbstractVector{<:Integer},
    magmom,
    spintype::SpinType,
    symprec::Real,
)
    if spintype isa Unpolarized
        # No magnetic constraints - both forward and backward are compatible
        return (forward = true, backward = true)
    elseif spintype isa Collinear
        # Collinear: check if moments match after permutation
        # Allow for spin flip (time reversal)
        natoms = length(perm)
        forward_ok = true
        backward_ok = true

        for i = 1:natoms
            j = perm[i]
            m_i = magmom[i]
            m_j = magmom[j]

            # Forward: moments should match (|m_i - m_j| < symprec)
            if abs(m_i - m_j) >= symprec
                forward_ok = false
            end

            # Backward (time reversal): moments can flip (|m_i + m_j| < symprec)
            if abs(m_i + m_j) >= symprec
                backward_ok = false
            end
        end

        return (forward = forward_ok, backward = backward_ok)
    else  # NonCollinear
        # Non-collinear: apply rotation to magnetic moment vectors
        natoms = length(perm)
        forward_ok = true
        backward_ok = true

        for i = 1:natoms
            j = perm[i]
            m_i = magmom[i]
            m_j = magmom[j]

            # Forward: R * m_i should equal m_j
            rotated = rotation * m_i
            if norm(rotated - m_j) >= symprec
                forward_ok = false
            end

            # Backward (time reversal): R * m_i should equal -m_j
            if norm(rotated + m_j) >= symprec
                backward_ok = false
            end
        end

        return (forward = forward_ok, backward = backward_ok)
    end
end

function my_compute_radius(lattvec::AbstractMatrix, nrotations::Integer, nkpt_target::Integer)
    vol = abs(det(lattvec))
    npoints = nkpt_target * nrotations
    return (3.0 / (4π) * npoints * vol)^(1 / 3)
end


function my_get_spacegroup_info(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    @info "\nPass here 21 in my_get_spacegroup_info: "
    Natoms = size(positions, 1)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:Natoms]
    return my_get_spacegroup_info(lattvec, pos_vec, types; symprec)
end

function my_get_spacegroup_info(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    @info "\nPass here 33 in my_get_spacegroup_info"
    cell = Spglib.Cell(lattvec, positions, types)
    dataset = Spglib.get_dataset(cell, symprec)
    return (
        spacegroup_number = dataset.spacegroup_number,
        international_symbol = String(dataset.international_symbol),
    )
end


function my_compute_bounds(lattvec::AbstractMatrix, radius::Real)
    # Metric tensor G = L^T L
    metric = lattvec' * lattvec

    # Inverse metric G^{-1}
    invmetric = inv(metric)

    # bounds[k] = ceil(r × √(G⁻¹ₖₖ))
    bounds = Vector{Int}(undef, 3)
    for k = 1:3
        bounds[k] = ceil(Int, radius * sqrt(invmetric[k, k]))
    end

    return bounds
end


function my_lattice_points_in_sphere(lattvec, radius, bounds)
    metric = lattvec' * lattvec  # G = L^T L
    r2 = radius^2
    points = Vector{SVector{3,Int}}()
    sqnorms = Vector{Float64}()
    for I in CartesianIndices((
        (-bounds[1]):bounds[1],
        (-bounds[2]):bounds[2],
        (-bounds[3]):bounds[3],
    ))
        n = SVector{3,Int}(Tuple(I))
        norm_sq = n' * metric * n
        if norm_sq < r2  # Strict inequality, matching C++ code
            push!(points, n)
            push!(sqnorms, norm_sq)
        end
    end

    # Sort by squared norm
    perm = sortperm(sqnorms)
    return points[perm], sqnorms[perm]
end


function my_get_unique_rotations(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Determine magnetic type
    mtype = my_get_magmom_type(magmom)

    # Create Spglib cell
    cell = Spglib.Cell(lattvec, positions, types)

    # Get symmetry dataset
    dataset = Spglib.get_dataset(cell, symprec)

    # Get rotations and translations
    rotations = dataset.rotations
    translations = dataset.translations

    # Compute unique rotations with time-reversal
    unique_rots = Set{Matrix{Int}}()
    rot_to_crot = Dict{Matrix{Int},Matrix{Float64}}()

    # Prepare for Cartesian conversion
    # Correct formula: R_cart = L * R_frac * L^{-1}
    L = lattvec
    L_inv = inv(L)

    for (iop, rot) in enumerate(rotations)
        trans = translations[iop]

        # Get atomic permutation
        perm = my_find_permutation(lattvec, positions, types, rot, trans, symprec)

        # Convert rotation to Cartesian coordinates: R_cart = L * R * L^{-1}
        rot_cart = L * Float64.(rot) * L_inv

        # Check magnetic compatibility
        compat = my_determine_compatibility(rot_cart, perm, magmom, mtype, symprec)

        # Store rotation directly (NOT transposed)
        # The C++ code stores wrapper.transpose() where wrapper is already R^T
        # due to Eigen's column-major interpretation of row-major spglib data,
        # so C++ stores (R^T)^T = R. We should store R directly.
        rot_int = Matrix{Int}(rot)

        if compat.forward
            if rot_int ∉ unique_rots
                push!(unique_rots, rot_int)
                # Cartesian rotation: R_cart = L * R * L^{-1}
                rot_to_crot[rot_int] = rot_cart
            end
        end
        if compat.backward
            neg_rot = -rot_int
            if neg_rot ∉ unique_rots
                push!(unique_rots, neg_rot)
                # Cartesian rotation for backward (time reversal): -R_cart
                rot_to_crot[neg_rot] = -rot_cart
            end
        end
    end

    # Convert to vectors
    rots_vec = collect(unique_rots)
    crots_vec = [rot_to_crot[r] for r in rots_vec]

    return rots_vec, crots_vec
end

function my_get_unique_rotations(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    Natoms = size(positions, 1)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:Natoms]
    return my_get_unique_rotations(lattvec, pos_vec, types, magmom; symprec)
end


function my_compute_equivalence_classes(points, sqnorms, rotations, symprec = 1e-5)
    n = length(points)

    # Build a lookup table for points
    point_to_idx = Dict{SVector{3,Int},Int}()
    for (i, p) in enumerate(points)
        point_to_idx[p] = i
    end

    # Mapping: -1 means not yet assigned
    mapping = fill(-1, n)

    for i = 1:n
        point = points[i]
        sqnorm = sqnorms[i]

        # Points can only be equivalent if they have the same norm
        # Binary search for range of candidates with similar norm
        lo = (1.0 - symprec) * sqnorm
        hi = (1.0 + symprec) * sqnorm
        lbound = searchsortedfirst(sqnorms, lo)
        ubound = min(searchsortedlast(sqnorms, hi), i - 1)  # Only earlier points

        # Check each rotation
        for R in rotations
            image = SVector{3,Int}(R * point)
            j = get(point_to_idx, image, 0)
            if j > 0 && lbound <= j <= ubound
                # Found a match - find its representative
                k = j
                while mapping[k] != k
                    k = mapping[k]
                end
                mapping[i] = k
                break
            end
        end

        if mapping[i] != -1
            # Already found a class, skip remaining rotations
            continue
        end

        # No equivalence found - start a new class
        if mapping[i] == -1
            mapping[i] = i
        end
    end

    return mapping
end

function my_calc_sphere_quotient_set(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom,
    radius::Real,
    bounds::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    Natoms = size(positions, 1)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:Natoms]
    return my_calc_sphere_quotient_set(
        lattvec,
        pos_vec,
        types,
        magmom,
        radius,
        bounds;
        symprec,
    )
end

function my_calc_sphere_quotient_set(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom,
    radius::Real,
    bounds::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    # Get unique rotations (integer matrices for reciprocal space)
    rotations, _ = my_get_unique_rotations(lattvec, positions, types, magmom; symprec)

    # Compute equivalence classes
    points, sqnorms = my_lattice_points_in_sphere(lattvec, radius, bounds)
    mapping = my_compute_equivalence_classes(points, sqnorms, rotations, symprec)

    # Group by equivalence - return as Vector of Vectors (easier to compare)
    classes = Dict{Int,Vector{SVector{3,Int}}}()
    for (i, rep) in enumerate(mapping)
        if !haskey(classes, rep)
            classes[rep] = SVector{3,Int}[]
        end
        push!(classes[rep], points[i])
    end

    # Sort by representative index to get consistent ordering
    reps = sort(collect(keys(classes)))
    return [classes[r] for r in reps]
end



function my_get_equivalences(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom,
    nkpt_target::Integer;
    symprec::Real = 1e-5,
)

    @info "\nPass here 46 in my_get_equivalences"

    # Ensure positions are natoms×3
    if size(positions, 1) == 3 && size(positions, 2) != 3
        positions = positions'
    end

    # Convert to Vector of SVectors for internal functions
    Natoms = size(positions, 1)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:Natoms]

    # Get number of rotations
    nrot = my_calc_nrotations(lattvec, pos_vec, types, magmom; symprec)

    # Compute radius for target number of equivalences
    radius = my_compute_radius(lattvec, nrot, nkpt_target)

    # Compute bounds
    bounds = my_compute_bounds(lattvec, radius)

    # Compute equivalence classes
    equivalences =
        my_calc_sphere_quotient_set(lattvec, pos_vec, types, magmom, radius, bounds; symprec)

    return equivalences, radius, nrot
end


function my_calc_nrotations(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Convert positions to vector of SVectors for Spglib
    Natoms = size(positions, 1)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:Natoms]
    return my_calc_nrotations(lattvec, pos_vec, types, magmom; symprec)
end

function my_calc_nrotations(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Determine magnetic type
    mtype = my_get_magmom_type(magmom)

    # Create Spglib cell
    cell = Spglib.Cell(lattvec, positions, types)

    # Get symmetry dataset
    dataset = Spglib.get_dataset(cell, symprec)

    # Get rotations and translations
    rotations = dataset.rotations
    translations = dataset.translations

    # Compute unique rotations with time-reversal
    unique_rots = Set{Matrix{Int}}()

    # Prepare for Cartesian conversion
    # Correct formula: R_cart = L * R_frac * L^{-1}
    L = lattvec
    L_inv = inv(L)

    for (iop, rot) in enumerate(rotations)
        trans = translations[iop]

        # Get atomic permutation
        perm = my_find_permutation(lattvec, positions, types, rot, trans, symprec)

        # Convert rotation to Cartesian coordinates
        rot_cart = L * Float64.(rot) * L_inv

        # Check magnetic compatibility
        compat = my_determine_compatibility(rot_cart, perm, magmom, mtype, symprec)

        # Store rotation directly (NOT transposed)
        # The C++ code stores wrapper.transpose() where wrapper is already R^T
        # due to Eigen's column-major interpretation of row-major spglib data,
        # so C++ stores (R^T)^T = R. We should store R directly.
        rot_int = Matrix{Int}(rot)

        if compat.forward
            push!(unique_rots, rot_int)
        end
        if compat.backward
            # Time reversal: -R
            push!(unique_rots, -rot_int)
        end
    end

    return length(unique_rots)
end


function my_find_permutation(
    lattvec::AbstractMatrix,
    positions::AbstractVector,
    types::AbstractVector{<:Integer},
    rotation::AbstractMatrix,
    translation::AbstractVector,
    symprec::Real,
)
    natoms = length(positions)
    perm = Vector{Int}(undef, natoms)

    for i = 1:natoms
        # Apply symmetry operation: R * pos + t
        pos_i = positions[i]
        new_pos = rotation * pos_i + translation

        # Wrap to [0, 1) range
        new_pos = mod.(new_pos, 1.0)

        # Find matching atom
        found = false
        for j = 1:natoms
            if types[i] != types[j]
                continue
            end

            # Compute difference in fractional coordinates
            diff = new_pos - positions[j]
            # Account for periodic boundary
            diff = mod.(diff .+ 0.5, 1.0) .- 0.5

            # Convert to Cartesian and check distance
            cart_diff = lattvec * diff
            if norm(cart_diff) < symprec
                perm[i] = j
                found = true
                break
            end
        end

        if !found
            error("Could not find permutation for atom $i")
        end
    end

    return perm
end


struct MyFourierInterpolator <: AbstractInterpolator
    coeffs::Matrix{ComplexF64}
    equivalences::Vector{<:AbstractMatrix{<:Integer}}
    lattvec::SMatrix{3,3,Float64}
end

function MyFourierInterpolator(kpoints, energies, equivalences, lattvec)
    coeffs = my_fitde3D(kpoints, energies, equivalences, lattvec)
    lattvec_static = SMatrix{3,3,Float64}(lattvec)
    return MyFourierInterpolator(coeffs, equivalences, lattvec_static)
end

function my_compute_phase_factors(kpoints, equivalences)
    nk = size(kpoints, 1)
    neq = length(equivalences)
    phase = zeros(ComplexF64, nk, neq)
    tpii = 2π * im
    for (j, equiv) in enumerate(equivalences)
        nstar = size(equiv, 1)
        # Vectorized: matrix multiply (nk, 3) * (3, nstar) -> (nk, nstar)
        dot_products = kpoints * equiv'
        # exp and sum over stars (dims=2)
        phase_sum = sum(exp.(tpii .* dot_products), dims = 2)
        phase[:, j] = vec(phase_sum) ./ nstar
    end
    return phase
end

function my_compute_regularization_weights(equivalences, lattvec; C1 = 0.75, C2 = 0.75)
    # Get representative vectors and their norms
    Rvecs = [equiv[1, :] for equiv in equivalences]
    norms = [norm(lattvec' * R) for R in Rvecs]
    # Normalize by first non-zero norm
    norm1 = norms[2]  # First non-zero (index 2, since index 1 is origin)
    X2 = (norms ./ norm1) .^ 2
    rhoi = @. 1.0 / ((1.0 - C1 * X2)^2 + C2 * X2^3)
    rhoi[1] = 0.0  # Origin has zero weight
    return rhoi
end


function my_fitde3D(kpoints, energies, equivalences, lattvec)

    @info "Pass here 521 in my_fitde3d"

    nk = size(kpoints, 1)
    nbands = size(energies, 1)
    neq = length(equivalences)

    # Compute phase factors and regularization weights
    phase = my_compute_phase_factors(kpoints, equivalences)
    rhoi = my_compute_regularization_weights(equivalences, lattvec)

    # Energy differences (relative to last k-point)
    # De has shape (nbands, nk-1) in Julia, corresponds to De.T in Python
    De = energies[:, 1:(end-1)] .- energies[:, end:end]

    # Phase differences: phaseR (nk-1, neq)
    phaseR = phase[1:(end-1), :] .- phase[end:end, :]

    # Build regularized matrix (matches Python exactly)
    # Python: Hmat = (phaseR[:, 1:] @ (phaseR[:, 1:] * rhoi[1:]).conj().T).real
    # phaseR[:, 2:end] is (nk-1, neq-1)
    # Element-wise multiply with rhoi (broadcast along rows)
    weighted_phaseR = phaseR[:, 2:end] .* rhoi[2:end]'  # (nk-1, neq-1)

    # Hmat = phaseR[:, 2:end] @ weighted_phaseR.conj().T = (nk-1, nk-1)
    Hmat = real(phaseR[:, 2:end] * conj(weighted_phaseR)')

    # Solve least squares: Hmat @ rlambda = De.T
    # Python: rlambda = sp.linalg.lstsq(Hmat, De)[0]  (SVD, rank-deficient OK)
    # De.T in Python is (nk-1, nbands), our De is (nbands, nk-1)
    # So we need De' for the solve.
    # Julia's `\` on a square matrix uses an LU factorization, which raises
    # SingularException for a rank-deficient Hmat (high-symmetry small cells
    # where the source k-mesh is sparse relative to the equivalence stars).
    # `pinv` (Moore-Penrose, SVD-based) reproduces scipy.linalg.lstsq's
    # minimum-norm least-squares solution and is identical to `\` for the
    # full-rank case up to round-off.
    rlambda = pinv(Hmat) * De'  # (nk-1, nbands)

    # Recover coefficients
    # Python: coeffs = rhoi * (rlambda.T @ phaseR)
    # rlambda.T is (nbands, nk-1), phaseR is (nk-1, neq)
    # rlambda.T @ phaseR = (nbands, neq)
    coeffs = rhoi' .* (rlambda' * phaseR)  # (nbands, neq)

    # First coefficient: E₀ = E_ref - Σ c_R exp(2πi k_ref · R)
    # Python: coeffs[:, 0] = ene.T[-1] - coeffs[:, 1:] @ phase[-1, 1:]
    coeffs[:, 1] = energies[:, end] - coeffs[:, 2:end] * phase[end, 2:end]

    return coeffs
end


struct MyInterpolationResult
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
function MyInterpolationResult(
    interp::MyFourierInterpolator;
    atoms = nothing,
    metadata = Dict{String,Any}(),
)
    # Convert equivalences to plain Matrix{Int} for serialization
    equivs = [Matrix{Int}(eq) for eq in interp.equivalences]
    lattvec = Matrix{Float64}(interp.lattvec)

    meta = Dict{String,Any}(
        "version" => "0.1.0",
        "created" => string(Dates.now()),
        "nbands" => size(interp.coeffs, 1),
        "neq" => length(interp.equivalences),
        "spintype" => "Unpolarized",  # v0.2 forward compatibility
    )
    merge!(meta, metadata)

    return MyInterpolationResult(interp.coeffs, equivs, lattvec, atoms, meta)
end



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
    sginfo = my_get_spacegroup_info(lattvec, positions', types; symprec)
    verbose && println(
        "  Space group: $(sginfo.spacegroup_number) ($(sginfo.international_symbol))",
    )

    # 4. Compute equivalences
    verbose && println("Computing equivalences for ~$nkpt_target k-points...")
    equivalences, radius, nrot =
        my_get_equivalences(lattvec, positions', types, magmom, nkpt_target; symprec)
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
    interp = MyFourierInterpolator(kpts, ebands, equiv_matrices, lattvec)

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

    result = MyInterpolationResult(interp; atoms = atoms, metadata = metadata)
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
    serialize("TEMP_interp.jldat", interp)
    return
end

