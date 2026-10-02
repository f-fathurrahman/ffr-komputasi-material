# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/bandlib.py

using LinearAlgebra
using Base.Threads: @threads, nthreads, threadid

# Fermi-Dirac distribution thresholds
const FD_XMAX = 18.4206807340       # fFD(x) = 1e-8 threshold for numerical stability
const FD_XMAX_GAP = 6.9067547786    # fFD(x) = 1e-3 threshold for gap detection

#=
    fermi_dirac(ε::Real, μ::Real, kT::Real) -> Float64

Compute Fermi-Dirac distribution f(ε) = 1/(exp((ε-μ)/kT) + 1) for a single energy.
=#
function fermi_dirac(ε::Real, μ::Real, kT::Real)
    x = (ε - μ) / kT
    if x < -FD_XMAX
        return 1.0
    elseif x > FD_XMAX
        return 0.0
    else
        return 1.0 / (exp(x) + 1.0)
    end
end

#=
    fermi_dirac!(f, ε, μ, kT)

Compute Fermi-Dirac distribution f(ε) = 1/(exp((ε-μ)/kT) + 1) in-place.
=#
function fermi_dirac!(f::AbstractArray, ε::AbstractArray, μ::Real, kT::Real)
    @inbounds for i in eachindex(ε, f)
        f[i] = fermi_dirac(ε[i], μ, kT)
    end
    return f
end

#=
    fermi_dirac(ε::AbstractArray, μ, kT)

Compute Fermi-Dirac distribution f(ε) = 1/(exp((ε-μ)/kT) + 1) for an array.
=#
function fermi_dirac(ε::AbstractArray, μ::Real, kT::Real)
    f = similar(ε, Float64)
    return fermi_dirac!(f, ε, μ, kT)
end

#=
    dfermi_dirac_de(ε::Real, μ::Real, kT::Real) -> Float64

Compute derivative of Fermi-Dirac distribution ∂f/∂ε for a single energy.
=#
function dfermi_dirac_de(ε::Real, μ::Real, kT::Real)
    x = (ε - μ) / kT
    if abs(x) > FD_XMAX
        return 0.0
    else
        c = cosh(0.5 * x)
        return -0.25 / (c * c * kT)
    end
end

#=
    dfermi_dirac_de!(df, ε, μ, kT)

Compute derivative of Fermi-Dirac distribution ∂f/∂ε in-place.
=#
function dfermi_dirac_de!(df::AbstractArray, ε::AbstractArray, μ::Real, kT::Real)
    @inbounds for i in eachindex(ε, df)
        df[i] = dfermi_dirac_de(ε[i], μ, kT)
    end
    return df
end

#=
    dfermi_dirac_de(ε::AbstractArray, μ, kT)

Compute derivative of Fermi-Dirac distribution ∂f/∂ε for an array.
=#
function dfermi_dirac_de(ε::AbstractArray, μ::Real, kT::Real)
    df = similar(ε, Float64)
    return dfermi_dirac_de!(df, ε, μ, kT)
end

#=
    compute_dos(eband, erange, npts; weights=nothing)

Compute density of states using histogram.

# Arguments
- `eband`: Band energies (nbands × nkpts)
- `erange`: (emin, emax) energy range
- `npts`: Number of energy bins
- `weights`: Optional weights for weighted DOS

# Returns
- `centers`: Energy bin centers
- `dos`: Density of states
=#
function compute_dos(eband::AbstractMatrix, erange::Tuple, npts::Int; weights = nothing)
    edges = range(erange[1], erange[2], length = npts+1)
    de = step(edges)
    centers = (collect(edges[1:(end-1)]) .+ collect(edges[2:end])) ./ 2

    if isnothing(weights)
        h = fit(Histogram, vec(eband), edges)
    else
        h = fit(Histogram, vec(eband), Weights(vec(weights)), edges)
    end

    nkpt = size(eband, 2)
    dos = h.weights ./ nkpt ./ de

    return centers, dos
end

#=
    BTPDOS(eband, vvband; erange=nothing, npts=nothing)

Compute DOS and transport DOS (velocity-weighted).

# Returns
- `epsilon`: Energy bins
- `dos`: Density of states
- `vvdos`: Transport DOS [3, 3, npts]
=#
function BTPDOS(
    eband::AbstractMatrix,
    vvband::AbstractArray;
    erange = nothing,
    npts = nothing,
)
    nbands, nkpts = size(eband)

    if isnothing(erange)
        # Match Python: no padding, just min/max
        erange = (minimum(eband), maximum(eband))
    end

    if isnothing(npts)
        npts = 1000
    end

    # Standard DOS
    epsilon, dos = compute_dos(eband, erange, npts)

    # Transport DOS (v⊗v weighted)
    vvdos = zeros(3, 3, npts)
    for i = 1:3, j = i:3
        weights = @view vvband[:, i, j, :]
        _, vvdos[i, j, :] = compute_dos(eband, erange, npts; weights = weights)
        if i != j
            vvdos[j, i, :] = vvdos[i, j, :]  # Symmetric
        end
    end

    return epsilon, dos, vvdos
end

#=
    fermi_integrals(epsilon, dos, sigma, μ_range, T_range; dosweight=2.0)

Compute Fermi integrals L0, L1, L2.

# Returns
- `N`: Electron count [nT, nμ]
- `L0`: σ(ε) integral [nT, nμ, 3, 3]
- `L1`: (ε-μ)σ(ε) integral [nT, nμ, 3, 3]
- `L2`: (ε-μ)²σ(ε) integral [nT, nμ, 3, 3]
=#
function fermi_integrals(epsilon, dos, sigma, μ_range, T_range; dosweight = 2.0)
    nT = length(T_range)
    nμ = length(μ_range)
    npts = length(epsilon)
    de = epsilon[2] - epsilon[1]

    N = zeros(nT, nμ)
    L0 = zeros(nT, nμ, 3, 3)
    L1 = zeros(nT, nμ, 3, 3)
    L2 = zeros(nT, nμ, 3, 3)

    # Pre-allocate temporary arrays per thread
    f_tls = [zeros(npts) for _ = 1:nthreads()]
    df_tls = [zeros(npts) for _ = 1:nthreads()]
    int0_tls = [zeros(npts) for _ = 1:nthreads()]
    intn_tls = [zeros(npts) for _ = 1:nthreads()]

    # Create list of (temperature, chemical potential) pairs
    param_space = collect(Iterators.product(enumerate(T_range), enumerate(μ_range)))

    # Parallelize loop
    @threads for ((iT, T), (iμ, μ)) in param_space
        tid = threadid()
        f = f_tls[tid]
        df = df_tls[tid]
        int0 = int0_tls[tid]
        intn = intn_tls[tid]

        kT = KB_AU * T  # T in Kelvin, KB_AU in Ha/K

        # Call in-place functions to avoid memory allocation
        fermi_dirac!(f, epsilon, μ, kT)
        dfermi_dirac_de!(df, epsilon, μ, kT)

        s = 0.0
        @simd for i in eachindex(dos, f)
            s += dos[i] * f[i]
        end
        N[iT, iμ] = -dosweight * s * de

        @. int0 = -dosweight * df
        for i = 1:3, j = i:3
            # Element-wise operations using @.
            @. intn = int0 * sigma[i, j, :]
            L0[iT, iμ, i, j] = sum(intn) * de

            @. intn *= (epsilon - μ)
            L1[iT, iμ, i, j] = -sum(intn) * de

            @. intn *= (epsilon - μ)
            L2[iT, iμ, i, j] = sum(intn) * de
        end
    end

    # Apply symmetry
    for iT = 1:nT, iμ = 1:nμ
        for i = 1:3, j = (i+1):3
            L0[iT, iμ, j, i] = L0[iT, iμ, i, j]
            L1[iT, iμ, j, i] = L1[iT, iμ, i, j]
            L2[iT, iμ, j, i] = L2[iT, iμ, i, j]
        end
    end

    return N, L0, L1, L2
end

#=
    _safe_inverse(A::AbstractMatrix)

Compute inverse of symmetric positive-definite matrix with fallback.

Uses Cholesky decomposition for speed when matrix is well-conditioned,
falls back to Moore-Penrose pseudoinverse for singular/ill-conditioned cases.
=#
function _safe_inverse(A::AbstractMatrix)
    try
        return inv(cholesky(Symmetric(A)))
    catch
        return pinv(A)
    end
end

#=
    calc_onsager_coefficients(L0, L1, L2, T_range, vuc)

Calculate Onsager transport coefficients.

Unit conversions follow Python BoltzTraP2 conventions:
- σ is in S/(m·s) (conductivity/τ)
- S is in V/K (Seebeck coefficient)
- κ is in W/(m·K·s) (thermal conductivity/τ)

# Returns
- `σ`: Electrical conductivity [nT, nμ, 3, 3]
- `S`: Seebeck coefficient [nT, nμ, 3, 3]
- `κ`: Thermal conductivity [nT, nμ, 3, 3]
=#
function calc_onsager_coefficients(L0, L1, L2, T_range, vuc)
    nT, nμ = size(L0)[1:2]

    σ = similar(L0)
    S = similar(L0)
    κ = similar(L0)

    # Pre-compute T-independent factor (optimization: multiply instead of divide)
    inv_sigma_vuc = 1.0 / (SIGMA_CONV * vuc)

    for iT = 1:nT, iμ = 1:nμ
        T = T_range[iT]

        # Pre-compute T-dependent factors
        inv_seebeck_T_vuc = 1.0 / (T * SEEBECK_CONV * vuc)
        inv_kappa_T_vuc = 1.0 / (T * KAPPA_CONV * vuc)

        # L0 → σ/τ with BoltzTraP2-compatible unit conversion
        # Python: L11 = L0 / (Siemens / (Meter * Second)) / vuc
        L11 = @view(L0[iT, iμ, :, :]) * inv_sigma_vuc
        σ[iT, iμ, :, :] = L11

        # L1 → for Seebeck calculation
        # Python: L12 = L1 / T / (Volt * Siemens / (Meter * Second)) / vuc
        L12 = @view(L1[iT, iμ, :, :]) * inv_seebeck_T_vuc

        # Seebeck: S = σ⁻¹ L12
        # Use Cholesky when possible (faster), fallback to pinv for singular cases
        σ_inv = _safe_inverse(L11)
        S[iT, iμ, :, :] = σ_inv * L12

        # L2 → for thermal conductivity
        # Python: L22 = L2 / T / (Volt * Joule * Siemens / (Meter * Second * Coulomb)) / vuc
        L22 = @view(L2[iT, iμ, :, :]) * inv_kappa_T_vuc

        # Thermal conductivity: κ = L22 - T * σ * S * S
        # Since S = σ⁻¹ * L12, this becomes: κ = L22 - T * L12 * σ⁻¹ * L12
        κ[iT, iμ, :, :] = L22 - T * L12 * σ_inv * L12
    end

    return σ, S, κ
end

"""
    calc_cv(epsilon, dos, μ_range, T_range; dosweight=2.0) -> Matrix{Float64}

Compute the electronic contribution to the heat capacity.

# Arguments

  - `epsilon`: energy grid (npts,), in Hartree
  - `dos`: density of states (npts,), from `BTPDOS` or `TransportResult.dos`
  - `μ_range`: chemical potential values (nμ,), in Hartree
  - `T_range`: temperature values (nT,), in Kelvin
  - `dosweight`: spin degeneracy (2.0 for unpolarized, 1.0 for polarized)

# Returns

Matrix of shape (nT, nμ) with heat capacity in SI units (J/K).

# Example

```julia
# From TransportResult
transport = run_integrate(\"si_interp.jld2\"; temperatures = [300.0])
cv = calc_cv(transport.epsilon, transport.dos, transport.mu_values, transport.temperatures)
```
"""
function calc_cv(epsilon, dos, μ_range, T_range; dosweight = 2.0)
    nT = length(T_range)
    nμ = length(μ_range)
    cv = Matrix{Float64}(undef, nT, nμ)
    de = epsilon[2] - epsilon[1]

    for (iT, T) in enumerate(T_range)
        kBT = T * KB_AU  # KB_AU from units.jl
        for (iμ, μ) in enumerate(μ_range)
            df = dfermi_dirac_de(epsilon, μ, kBT)
            # cv = k_B_SI * (-dosweight * ∫ (ε-μ)²/kBT * (∂f/∂ε) * dos dε)
            # Following Python BoltzTraP2: -(dosweight * dos * dFDde * (ε-μ)² / kBT) * de
            integrand = -dosweight .* dos .* df .* (epsilon .- μ) .^ 2 ./ kBT
            cv[iT, iμ] = KB_SI * sum(integrand) * de
        end
    end

    return cv
end

#=
    calc_transport_coefficients(eband, vvband, T_range, μ_range, vuc)

High-level function to compute all transport coefficients.
=#
function calc_transport_coefficients(eband, vvband, T_range, μ_range, vuc)
    epsilon, dos, vvdos = BTPDOS(eband, vvband)
    N, L0, L1, L2 = fermi_integrals(epsilon, dos, vvdos, μ_range, T_range)
    σ, S, κ = calc_onsager_coefficients(L0, L1, L2, T_range, vuc)
    return (σ = σ, S = S, κ = κ, N = N, epsilon = epsilon, dos = dos)
end
