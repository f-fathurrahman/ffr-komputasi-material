using BoltzTraP
using BoltzTraP: FourierInterpolator, getBTPbands, BTPDOS
using BoltzTraP: KB_AU, BOHR_TO_ANG, HA_TO_EV
using BoltzTraP: solve_for_mu, fermi_dirac!, dfermi_dirac_de!, calc_onsager_coefficients
using StaticArrays
using LinearAlgebra: det
using Plots


function debug_interpolate()

    # Path to Si VASP data (directory containing vasprun.xml and POSCAR)
    datadir = "./data_Si_vasp"

    # Output file stem
    stem = "TEMP_DATA"

    # Step 1: Interpolate band structure
    println("Step 1: Interpolating band structure...")
    interp = run_interpolate(datadir; kpoints = 5000, verbose = true)
    save_interpolation(joinpath(@__DIR__, stem * "_interp.jld2"), interp)

    println("  Equivalence classes: $(length(interp.equivalences))")
    println("  Bands: $(size(interp.coeffs, 1))")
    
    return
end


function my_run_integrate(
    interp::InterpolationResult;
    temperatures::AbstractVector{<:Real} = [300.0],
    output::Union{String,Nothing} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)
    @debug "run_integrate called" temperatures output bins scissor verbose

    sys = TransportSystem(interp)
    fi = FourierInterpolator(
        interp.coeffs,
        interp.equivalences,
        SMatrix{3,3,Float64}(interp.lattvec),
    )

    result = my_solve_transport(UniformMesh(), fi, sys; temperatures, bins, scissor, verbose)

    # Restore source provenance from the interpolation metadata.
    result.metadata["source"] = get(interp.metadata, "source_file", "unknown")

    if !isnothing(output)
        verbose && println("Saving to $output...")
        save_integrate(output, result)
    end

    verbose && println("Done.")
    return result
end


function my_fermi_integrals(epsilon, dos, sigma, μ_range, T_range; dosweight = 2.0)
    nT = length(T_range)
    nμ = length(μ_range)
    npts = length(epsilon)
    de = epsilon[2] - epsilon[1]

    N = zeros(nT, nμ)
    L0 = zeros(nT, nμ, 3, 3)
    L1 = zeros(nT, nμ, 3, 3)
    L2 = zeros(nT, nμ, 3, 3)

    # Pre-allocate temporary arrays per thread
    f_tls = zeros(npts)
    df_tls = zeros(npts)
    int0_tls = zeros(npts)
    intn_tls = zeros(npts)

    # Create list of (temperature, chemical potential) pairs
    param_space = collect(Iterators.product(enumerate(T_range), enumerate(μ_range)))

    # Parallelize loop
    for ((iT, T), (iμ, μ)) in param_space
        #tid = threadid()
        # ffr: just assign to new names
        f = f_tls
        df = df_tls
        int0 = int0_tls
        intn = intn_tls

        kT = KB_AU * T  # T in Kelvin, KB_AU in Ha/K

        # Call in-place functions to avoid memory allocation
        fermi_dirac!(f, epsilon, μ, kT)
        dfermi_dirac_de!(df, epsilon, μ, kT)

        s = 0.0
        for i in eachindex(dos, f)
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


# Concrete dispatch: UniformMesh + FourierInterpolator
function my_solve_transport(
    sampling::UniformMesh,
    interp::FourierInterpolator,
    sys::TransportSystem;
    temperatures::AbstractVector{<:Real} = [300.0],
    mur::Union{Nothing,AbstractVector{<:Real}} = nothing,
    bins::Int = 0,
    scissor::Union{Nothing,Real} = nothing,
    verbose::Bool = false,
)::TransportResult
    nelect = sys.nelect
    dosweight = sys.dosweight
    fermi_dft = sys.fermi

    # Step 1: bands via FFT (Fourier-specific)
    eband, vvband = getBTPbands(interp.coeffs, interp.equivalences, interp.lattvec)
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
    N, L0, L1, L2 = my_fermi_integrals(epsilon, dos, vvdos, μ_range, Tr; dosweight)

    # Step 6: Onsager coefficients (volume in Å³ matches Python BoltzTraP2)
    lattvec_ang = sys.lattice * BOHR_TO_ANG
    vuc = abs(det(lattvec_ang))
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
        "vuc_ang3" => vuc,
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




function debug_transport()

    interp = load_interpolation("TEMP_DATA_interp.jld2")

    @debug "Interpolation data loaded" nbands=size(interp.coeffs, 1) nequiv=length(
        interp.equivalences,
    )
    # Step 2: Compute transport coefficients
    println("\nStep 2: Computing transport coefficients...")
    temperatures = [200.0, 300.0, 400.0, 500.0]
    transport = my_run_integrate(interp; temperatures = temperatures, verbose = true)

    stem = "TEMP_DATA"
    save_integrate(joinpath(@__DIR__, stem * "_transport.jld2"), transport)

    # Step 3: Print results
    println("\nResults:")
    println("  Temperatures: $(transport.temperatures) K")
    println("  Chemical potential points: $(length(transport.mu_values))")
    println("  Tensor shape (σ): $(size(transport.sigma))")

    # Step 4: Plot transport coefficients (S, σ, κ vs μ at T=300K)
    fermi_dft_eV = transport.metadata["fermi_dft_eV"]
    mu = transport.mu_values .- fermi_dft_eV
    iT = findfirst(t -> isapprox(t, 300.0, atol = 1.0), transport.temperatures)
    p1 = plot(
        mu,
        transport.seebeck[1, 1, iT, :] .* 1e6;
        title = "T = 300 K",
        ylabel = "S_xx (μV/K)",
        ylims = (-1100, 1100),
        xlabel = "",
        xformatter = _->"",
        legend = false,
        linewidth = 2,
    )
    p2 = plot(
        mu,
        transport.sigma[1, 1, iT, :];
        ylabel = "σ_xx/τ (S/m)",
        yscale = :log10,
        ylims = (1e14, 1e21),
        xlabel = "",
        xformatter = _->"",
        legend = false,
        linewidth = 2,
    )
    p3 = plot(
        mu,
        transport.kappa[1, 1, iT, :];
        ylabel = "κ_xx/τ (W/m/K)",
        yscale = :log10,
        ylims = (1e10, 1e16),
        xlabel = "μ - E_F (eV)",
        legend = false,
        linewidth = 2,
    )
    p = plot(
        p1,
        p2,
        p3;
        layout = (3, 1),
        size = (700, 900),
        xlims = (-0.5, 0.5),
        left_margin = 5Plots.mm,
    )
    output_png = joinpath(@__DIR__, stem * "_transport_300K.png")
    savefig(p, output_png)
    println("Saved plot to $output_png")
end
