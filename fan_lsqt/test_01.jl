function test_main()
    # (1) Prepare some parameters
    Nx = 50000;
    Ny = 2;
    W = 0.0;
    E_max = 3.1; # energy scaling factor
    M = 1000; # number of Chebyshev moments
    E = collect(range(-3, stop=3, length=601)); # energy points
    dt = ones(10); # time steps

    H, V = init_H_and_V(Nx, Ny, W);

    ϕ = init_state(Nx*Ny);

    DOS = calc_dos(M, E_max, E/E_max, H/E_max, ϕ);

    VAC, sigma_from_VAC = calc_vac(M, E_max, dt*E_max, E/E_max, H/E_max, V, ϕ, DOS);

    MSD, sigma_from_MSD = calc_msd(M, E_max, dt*E_max, E/E_max, H/E_max, V/E_max, ϕ, DOS);

    #dt_scaled = dt*E_max;
    #E_scaled = E/E_max;
    #H_scaled = H/E_max

    v_F = sqrt.(VAC[1,:]) # Fermi velocity
    g = 0.5*Ny * DOS .* v_F  # conductance
    g *= 2π # from e^2/hbar to e^2/h

end