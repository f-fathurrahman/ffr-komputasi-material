using Chairmarks
using Revise

includet("mod_lsqt_01.jl")

includet("run_small.jl")
includet("bench_H_and_V.jl")

function bench_mult_H()
    Nx = 50_000;
    Ny = 2;
    W = 0.0;
    E_max = 3.1; # energy scaling factor
    M = 1000; # number of Chebyshev moments
    E = collect(range(-3, stop=3, length=601)); # energy points
    dt = ones(10); # time steps
    @info "bench H mult"
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    
    #@info "Out-place"
    #@b H * ϕ
    Hv = similar(ϕ)
    @info "In-place"
    @b mul!(Hv, H, ϕ)
end
#bench_mult_H()

function bench_mult_V()
    Nx = 50_000;
    Ny = 2;
    W = 0.0;
    E_max = 3.1; # energy scaling factor
    M = 1000; # number of Chebyshev moments
    E = collect(range(-3, stop=3, length=601)); # energy points
    dt = ones(10); # time steps
    @info "bench V mult"
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    @b V * ϕ
end
#bench_mult_V()

function bench_DOS()
    Nx = 50_000
    Ny = 2
    W = 0.0
    E_max = 3.1
    M = 1000
    E = collect(range(-3, stop=3, length=601))
    dt = ones(10); # time steps
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    E_scaled = E/E_max
    H_scaled = H/E_max
    @info "bench DOS"
    @b calc_dos(M, E_max, E_scaled, H_scaled, ϕ)
end
#bench_DOS()


function bench_moments()
    Nx = 50_000
    Ny = 2
    W = 0.0
    E_max = 3.1
    M = 1000
    E = collect(range(-3, stop=3, length=601))
    dt = ones(10); # time steps
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    H_scaled = H/E_max
    C = zeros(Float64, M)
    @info "bench moments in-place"
    @b calc_moments!(M, H_scaled, ϕ, ϕ, C)
end


function bench_DOS_fused()
    Nx = 50_000
    Ny = 2
    W = 0.0
    E_max = 3.1
    M = 1000
    E = collect(range(-3, stop=3, length=601))
    dt = ones(10); # time steps
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    E_scaled = E/E_max
    H_scaled = H/E_max
    @info "bench DOS fused"
    @b calc_dos_fused(M, E_max, E_scaled, H_scaled, ϕ)
end
#bench_DOS_fused()



function bench_VAC()
    Nx = 50_000
    Ny = 2
    W = 0.0
    E_max = 3.1
    M = 1000
    E = collect(range(-3, stop=3, length=601))
    dt = ones(10); # time steps
    
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    
    E_scaled = E/E_max;
    dt_scaled = dt*E_max;
    V_scaled = V/E_max;
    H_scaled = H/E_max;

    DOS = calc_dos(M, E_max, E_scaled, H_scaled, ϕ)
    @info "bench VAC"
    @b calc_vac(M, E_max, dt_scaled, E_scaled, H_scaled, V, ϕ, DOS);
end
#bench_VAC()


function bench_MSD()
    Nx = 50_000
    Ny = 2
    W = 0.0
    E_max = 3.1
    M = 1000
    E = collect(range(-3, stop=3, length=601))
    dt = ones(10); # time steps
    
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    
    E_scaled = E/E_max;
    dt_scaled = dt*E_max;
    V_scaled = V/E_max;
    H_scaled = H/E_max;

    DOS = calc_dos(M, E_max, E_scaled, H_scaled, ϕ)
    @info "bench MSD"
    @b calc_msd(M, E_max, dt_scaled, E_scaled, H_scaled, V_scaled, ϕ, DOS);
end
#bench_MSD()




function bench_evolve()
    Nx = 50_000
    Ny = 2
    W = 0.0
    E_max = 3.1
    M = 1000
    E = collect(range(-3, stop=3, length=601))
    dt = ones(10); # time steps
    
    H, V = init_H_and_V(Nx, Ny, W)
    ϕ = init_state(Nx*Ny)
    
    E_scaled = E/E_max;
    dt_scaled = dt*E_max;
    V_scaled = V/E_max;
    H_scaled = H/E_max;

    @info "bench evolve"
    @b evolve(H_scaled, dt_scaled[1], -1, ϕ)
end
#bench_evolve()

