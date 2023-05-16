using LinearAlgebra
using Printf

using Plots
using PlotThemes
theme(:dark)

include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_11pt.jl")

function build_Ham(N::Int64, Δx::Float64)
    D2 = build_D2_matrix_11pt(N, Δx)
    return -0.5*D2
end

# Using Ha unit: ħ = 1, m = 1


#function gaussian_wavepacket(x, t, σ, k0; x0=0.0)
#    Ω = σ^2
#    v = k0
#    oneo = (1 + im*Ω*t)
#    pw_f = exp( im*k0*(x - v*t - x0) ) # shift by x0 ?
#    gauss_f = exp( -σ^2/2 * (x - v*t - x0)^2 / oneo )
#    pre_f = sqrt(σ/( sqrt(π)*oneo) )
#    f = pre_f * gauss_f * pw_f
#    return f
#end


function gaussian_wavepacket(x, t, a0, p0)
    dc = 1.0 + 2.0im*t/a0^2
    pw_f = exp( im*p0*(x - 0.5*p0^2*t) )
    gauss_f = exp( -(x - p0*t)^2/(a0^2*dc) )
    pre_f = (2.0/(pi*a0^2 * dc^2))^0.25
    f = pre_f * gauss_f * pw_f
    return f
end


function main()

    a = -25.0
    b =  25.0

    σ = sqrt(1.0)
    k0 = 20.0
    x0 = 0.0

    k_grid_min = 2*k0 * 10
    Δx = 2π/k_grid_min
    println("Δx min = ", Δx)
    Nx = floor(Int64, (b - a)/Δx) + 2
    println("Calculated Nx = ", Nx)

    if Nx < 100
        Nx = 200
        println("Nx is set to ", Nx)
    end

    # NOTE: This is different from the one used in Varga
    x = zeros(Float64, Nx)
    Δx = (b - a)/(Nx - 1)
    for i in 1:Nx
        x[i] = a + (i-1)*Δx # calculate x grid point here
    end

    t_last = 0.5
    Nt = 2000
    Δt = t_last/(Nt - 1)
    
    println("Δt = ", Δt)

    println("k0*t_last = ", k0*t_last)

    H = build_Ham(Nx, Δx)

    U = inv(I + im*Δt/2 * H)*(I - im*Δt/2 * H)
    
    #U_RHS = I - im*Δt/2 * H
    #U_LHS = I + im*Δt/2 * H

    println("2π/k0 = ", 2π/k0)
    println("Δx    = ", Δx)

    println()
    println("k0     = ", k0)
    println("k grid = ", 2π/Δx)

    println("CFL = ", k0*Δt/Δx)

    # initial wave function
    t = 0.0
    ψ = zeros(ComplexF64, Nx)
    for i in 1:Nx
        ψ[i] = gaussian_wavepacket(x[i], t, σ, k0)
    end

    # Check norm
    @printf("Initial norm of the wavepacket: %18.10f\n", dot(ψ, ψ)*Δx)


    plot(x, abs.(ψ), label="Abs", linewidth=2, dpi=150)
    title!("Initial state")
    ylims!(0.0, 0.8)
    xlims!(a, b)
    savefig("IMG_initial.png")

    plot(x, real.(ψ), label="Re", linewidth=2, dpi=150)
    title!("Initial state")
    #ylims!(0.0, 0.8)
    xlims!(0.1*a, 0.1*b)
    savefig("IMG_initial_Re.png")

    t = 0.0
    ψ_ana = zeros(ComplexF64, Nx)
    ψ_new = zeros(ComplexF64, Nx)
    for itime in 1:Nt
        t += Δt
        @views ψ_new[:] = U * ψ[:]
        for i in 1:Nx
            ψ_ana[i] = gaussian_wavepacket(x[i], t, σ, k0)
        end
        ψ[:] = ψ_new[:]
    end
    println("Final t = ", t)

    plot(x, abs.(ψ_new), label="Numerical", linewidth=2, dpi=150)
    title!("Final state")
    ylims!(0.0, 0.8)
    xlims!(a, b)
    plot!(x, abs.(ψ_ana), label="Analytical", linewidth=2, dpi=150)
    savefig("IMG_final.png")

    @printf("After numerical propagation: norm of the wavepacket: %18.10f\n", dot(ψ,ψ)*Δx)

end
main()



#=
# Using analytical expression
t_last = Nt*Δt
nrm = 0.0
ψ_ana = zeros(ComplexF64, Nx)
for i in 1:Nx
    ψ_ana[i] = gaussian_wavepacket(x[i], t_last, σ, k0, x0=x0)
    nrm = nrm + real(conj(ψ_ana[i])*ψ_ana[i]*dx)
end

plot(x, abs.(ψ_ana), linewidth=2, dpi=150)
ylims!(0.0, 0.8)
savefig("IMG_itime_" * string(Nt) * ".png")

@printf("Analytical propagation: norm of the wavepacket: %18.10f\n", nrm)
=#
