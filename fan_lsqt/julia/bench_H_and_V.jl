function bench_H_and_V()
    Nx = 50_000;
    Ny = 2;
    W = 0.0;
    E_max = 3.1; # energy scaling factor
    M = 1000; # number of Chebyshev moments
    E = collect(range(-3, stop=3, length=601)); # energy points
    dt = ones(10); # time steps
    @info "bench H and V"
    @b init_H_and_V(Nx, Ny, W)
end
#bench_H_and_V()