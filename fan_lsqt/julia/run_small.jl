function run_small()
    @info "Run with small parameters"
    Nx = 10;
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
    @info "End of run with small parameters"
end
run_small()