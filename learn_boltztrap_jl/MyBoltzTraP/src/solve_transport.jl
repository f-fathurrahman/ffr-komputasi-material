# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

"""
    solve_transport(sampling, interp, sys; kwargs...) -> TransportResult

Method-dispatched entry point for transport calculations on the
`AbstractBZSampling` × `AbstractInterpolator` × `TransportSystem` axes.

Method overloads:

  - `::Any` fallback: throws `ArgumentError` naming the offending argument type
  - `::AbstractBZSampling, ::AbstractInterpolator, ::TransportSystem`:
    throws `ErrorException` for a sampling subtype that has no concrete method
  - concrete pair (e.g. `UniformMesh + FourierInterpolator`): actual computation

The concrete `UniformMesh + FourierInterpolator` method below mirrors the
existing `run_integrate(::InterpolationResult)` flow so its output is
bit-equal under default kwargs.
"""
function solve_transport(sampling, interp, sys; kwargs...)
    sampling isa AbstractBZSampling || throw(
        ArgumentError(
            "solve_transport: sampling must be ::AbstractBZSampling, got $(typeof(sampling)).",
        ),
    )
    interp isa AbstractInterpolator || throw(
        ArgumentError(
            "solve_transport: interp must be ::AbstractInterpolator, got $(typeof(interp)).",
        ),
    )
    sys isa TransportSystem || throw(
        ArgumentError(
            "solve_transport: sys must be ::TransportSystem, got $(typeof(sys)).",
        ),
    )
    throw(
        ArgumentError(
            "solve_transport: no concrete method for " *
            "($(typeof(sampling)), $(typeof(interp)), $(typeof(sys))).",
        ),
    )
end

# Subtype-without-impl fallback for AbstractBZSampling subtypes that do not
# provide a concrete method.
function solve_transport(
    sampling::AbstractBZSampling,
    interp::AbstractInterpolator,
    sys::TransportSystem;
    kwargs...,
)
    error(
        "solve_transport: not implemented for sampling type $(typeof(sampling)). " *
        "Concrete subtypes of AbstractBZSampling must define their own method.",
    )
end

# Concrete dispatch: UniformMesh + FourierInterpolator
function solve_transport(
    sampling::UniformMesh,
    interp::FourierInterpolator,
    sys::TransportSystem;
    temperatures::AbstractVector{<:Real} = [300.0],
    mur::Union{Nothing,AbstractVector{<:Real}} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)::TransportResult
    isnothing(sampling.nk) || @warn(
        "solve_transport: the Fourier path takes its integration grid from the " *
        "interpolation basis, so the requested UniformMesh size is ignored.",
        requested = sampling.nk,
        used = determine_fft_dims(interp.equivalences),
    )
    eband, vvband = getBTPbands(interp.coeffs, interp.equivalences, interp.lattvec)
    return _solve_transport_from_bands(
        eband,
        vvband,
        sys,
        sampling,
        interp;
        temperatures,
        mur,
        bins,
        scissor,
        verbose,
    )
end

# Concrete dispatch: UniformMesh + WannierInterpolator
function solve_transport(
    sampling::UniformMesh,
    interp::WannierInterpolator,
    sys::TransportSystem;
    temperatures::AbstractVector{<:Real} = [300.0],
    mur::Union{Nothing,AbstractVector{<:Real}} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)::TransportResult
    eband, vvband = getBTPbands_wannier(interp, sampling)
    return _solve_transport_from_bands(
        eband,
        vvband,
        sys,
        sampling,
        interp;
        temperatures,
        mur,
        bins,
        scissor,
        verbose,
    )
end

# Shared transport pipeline: DOS, μ-resolution, Fermi integrals, Onsager
# coefficients. Both Fourier and Wannier paths feed `(eband, vvband)` of
# matching shape — `(nbands, npts)` and `(nbands, 3, 3, npts)` — into
# this helper, so output is bit-equal to the prior monolithic Fourier
# implementation under default kwargs.
#=
    integration_grid(interp, sampling) -> Union{Nothing,NTuple{3,Int}}

The k-grid that the integration actually used. The Fourier path takes it from
the interpolation basis, whose star functions fix the smallest grid that
represents them exactly; the Wannier path takes it from the sampling, which
must carry an explicit `nk`.
=#
integration_grid(::AbstractInterpolator, sampling::AbstractBZSampling) = sampling.nk

integration_grid(interp::FourierInterpolator, ::UniformMesh) =
    determine_fft_dims(interp.equivalences)

function _solve_transport_from_bands(
    eband::AbstractMatrix{<:Real},
    vvband::AbstractArray{<:Real,4},
    sys::TransportSystem,
    sampling::AbstractBZSampling,
    interp::AbstractInterpolator;
    temperatures::AbstractVector{<:Real},
    mur::Union{Nothing,AbstractVector{<:Real}},
    bins::Int,
    scissor::Union{Nothing,Real},
    verbose::Bool,
)::TransportResult
    nelect = sys.nelect
    dosweight = sys.dosweight
    fermi_dft = sys.fermi

    nbands, npts = size(eband)

    # Step 1.5: optional scissor correction
    scissor_Ha = nothing
    if !isnothing(scissor)
        scissor_Ha = scissor * EV_TO_HA
        npts_dos_temp = bins > 0 ? bins : 500
        epsilon_temp, dos_temp, _ = BTPDOS(eband, vvband; npts = npts_dos_temp)
        eband = apply_scissor(epsilon_temp, dos_temp, nelect, eband, scissor_Ha; dosweight)
    end

    # Step 2: DOS and transport DOS
    npts_dos = bins > 0 ? bins : 500
    epsilon, dos, vvdos = BTPDOS(eband, vvband; npts = npts_dos)

    # Step 3: μ range — auto-generate from DOS grid, or use caller-supplied grid
    if isnothing(mur)
        margin = 9.0 * KB_AU * maximum(temperatures)
        μ_min = epsilon[1] + margin
        μ_max = epsilon[end] - margin
        if μ_min >= μ_max
            error("Energy window too narrow for requested temperatures")
        end
        μ_indices = findall(e -> e > μ_min && e < μ_max, epsilon)
        μ_range = epsilon[μ_indices]
    else
        μ_range = collect(Float64, mur)
    end

    # Step 4: μ for each T + refined Fermi (T=0)
    Tr = collect(Float64, temperatures)
    nT = length(Tr)
    μ0 = zeros(nT)
    for (iT, T) in enumerate(Tr)
        μ0[iT] = solve_for_mu(
            epsilon,
            dos,
            nelect,
            T;
            dosweight,
            refine = true,
            try_center = true,
        )
    end
    fermi_Ha =
        solve_for_mu(epsilon, dos, nelect, 0.0; dosweight, refine = true, try_center = true)

    # Step 5: Fermi integrals
    N, L0, L1, L2 = fermi_integrals(epsilon, dos, vvdos, μ_range, Tr; dosweight)

    # Step 6: Onsager coefficients (volume in Bohr³, atomic-units convention)
    vuc = abs(det(sys.lattice))
    σ, S, κ = calc_onsager_coefficients(L0, L1, L2, Tr, vuc)

    # Step 7: μ to eV for output
    μ_range_eV = μ_range .* HA_TO_EV

    # Step 8: result metadata + tensor reshape
    spintype_str = if sys.spintype isa Unpolarized
        "Unpolarized"
    elseif sys.spintype isa Collinear
        "Collinear"
    elseif sys.spintype isa NonCollinear
        "NonCollinear"
    else
        "Unknown"
    end

    result_metadata = Dict{String,Any}(
        "source" => "solve_transport",
        "sampling" => string(nameof(typeof(sampling))),
        "sampling_nk" => integration_grid(interp, sampling),
        "interpolator" => replace(string(nameof(typeof(interp))), "Interpolator" => ""),
        "nelect" => nelect,
        "dosweight" => dosweight,
        "spintype" => spintype_str,
        "fermi_Ha" => fermi_Ha,
        "fermi_eV" => fermi_Ha * HA_TO_EV,
        "fermi_dft_Ha" => fermi_dft,
        "fermi_dft_eV" => fermi_dft * HA_TO_EV,
        "mu0_Ha" => μ0,
        "mu0_eV" => μ0 .* HA_TO_EV,
        "vuc_bohr3" => vuc,
        "nbands" => nbands,
        "npts_fft" => npts,
        "npts_dos" => npts_dos,
    )
    if !isnothing(scissor)
        result_metadata["scissor_eV"] = scissor
        result_metadata["scissor_Ha"] = scissor_Ha
    end

    σ_out = permutedims(σ, (3, 4, 1, 2))
    S_out = permutedims(S, (3, 4, 1, 2))
    κ_out = permutedims(κ, (3, 4, 1, 2))

    dos_info = Dict{String,Any}(
        "epsilon_Ha" => epsilon,
        "epsilon_eV" => epsilon .* HA_TO_EV,
        "dos" => dos,
    )

    return TransportResult(Tr, μ_range_eV, σ_out, S_out, κ_out, dos_info, result_metadata)
end
