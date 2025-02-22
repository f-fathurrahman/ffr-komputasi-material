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

    C = find_moments(M, H/E_max, ϕ, ϕ);
    DOS = chebyshev_summation(M, C, E/E_max, E_max);
    #DOS = find_dos(M, E_max, E/E_max, H/E_max, ϕ);

end