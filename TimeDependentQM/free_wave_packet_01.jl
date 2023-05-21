using LinearAlgebra
using Printf

# Using Ha unit: ħ = 1, m = 1

function gaussian_wavepacket(x, t, σ, k0; x0=0.0)
    Ω = σ^2
    v = k0
    oneo = (1 + im*Ω*t)
    pw_f = exp( im*k0*(x - v*t) )
    gauss_f = exp( -σ^2/2 * (x - v*t - x0)^2 / oneo )
    pre_f = sqrt(σ/( sqrt(π)*oneo) )
    f = pre_f * gauss_f * pw_f
    return f
end


Nx = 200
Nt = 500

a = -25.0
b =  25.0

x = zeros(Float64, Nx)
dx = (b-a)/(Nx + 1)
for i in 1:Nx
    x[i] = a + i*dx # calculate x grid point here
end

σ = sqrt(1.0)
k0 = 20.0
x0 = -10.0
# initial wave function
t = 0.0
ψ = zeros(ComplexF64, Nx)
for i in 1:Nx
    ψ[i] = gaussian_wavepacket(x[i], t, σ, k0, x0=x0)
end

# Check norm
nrm = 0.0
for i in 1:Nx
    nrm = nrm + real(conj(ψ[i])*ψ[i]*dx)
end
@printf("Initial norm of the wavepacket: %18.10f\n", nrm)


using Plots
using PlotThemes
theme(:dark)

plot(x, abs.(ψ), linewidth=2, dpi=150)
ylims!(0.0, 0.8)
savefig("IMG_itime_" * string(0) * ".png")

# Using analytical expression
t_last = Nt*Δt
nrm = 0.0
for i in 1:Nx
    ψ[i] = gaussian_wavepacket(x[i], t_last, σ, k0, x0=x0)
    nrm = nrm + real(conj(ψ[i])*ψ[i]*dx)
end

plot(x, abs.(ψ), linewidth=2, dpi=150)
ylims!(0.0, 0.8)
savefig("IMG_itime_" * string(Nt) * ".png")

@printf("Analytical propagation: norm of the wavepacket: %18.10f\n", nrm)

