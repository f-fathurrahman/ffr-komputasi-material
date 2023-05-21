using LinearAlgebra
using Printf

using Plots
using PlotThemes
theme(:dark)

function gaussian_wavepacket(x, t, a0, p0)
    dc = 1.0 + 2.0im*t/a0^2
    pw_f = exp( im*p0*(x - 0.5*p0^2*t) )
    gauss_f = exp( -(x - p0*t)^2/(a0^2*dc) )
    pre_f = (2.0/(pi*a0^2 * dc^2))^0.25
    f = pre_f * gauss_f * pw_f
    return f
end


function main()
    Nx = 200
    Nt = 500

    h2m = 0.5
    xnu = 1.0

    a = -25.0
    b =  25.0


    x = zeros(Float64, Nx)
    Δx = (b - a)/(Nx + 1)
    u = zeros(Float64, Nx, Nx)
    for i in 1:Nx
        u[i,i] = 1.0
        x[i] = a + i*Δx # calculate x grid point here
    end

    H = zeros(Float64, Nx, Nx)
    t = h2m/Δx^2

    for i in 1:Nx
        
        H[i,i] = 205.0/72.0*t
        
        if i > 1
            H[i,i-1] = -8.0/5.0*t
            H[i-1,i] = -8.0/5.0*t
        end
    
        if i > 2
            H[i,i-2] = 1.0/5.0*t
            H[i-2,i] = 1.0/5.0*t
        end
    
        if i > 3
            H[i,i-3] = -8.0/315.0*t
            H[i-3,i] = -8.0/315.0*t
        end
    
        if i > 4
            H[i,i-4] = 1.0/560.0*t
            H[i-4,i] = 1.0/560.0*t
        end
    end
    
    Hc = zeros(ComplexF64, Nx, Nx)
    Hc[:,:] = H[:,:]
   
    Δt = 0.001
    am = u - 0.50im*Hc*Δt
    zo = u + 0.50im*Hc*Δt

    bm = inv(zo)
    cm = am*bm


    a0 = 2.0 # must be larger or equal to 1 ?
    p0 = 10.0
    # initial wave function
    t = 0.0
    c = zeros(ComplexF64, Nx)
    for i in 1:Nx
        c[i] = gaussian_wavepacket(x[i], t, a0, p0)
    end

    # Check norm
    @printf("Initial norm of the wavepacket: %18.10f\n", dot(c,c)*Δx)

    # Numerical propagation
    plot(x, abs.(c), linewidth=2, dpi=150)
    savefig("IMG_orig_itime_" * string(0) * ".png")
    for itime in 1:Nt
        @views c[:] = cm * c[:]
    end

    plot(x, abs.(c), label="Numerical", linewidth=2, dpi=150)

    @printf("After numerical propagation: norm of the wavepacket: %18.10f\n", dot(c,c)*Δx)

    t_last = Nt*Δt
    ψ_ana = zeros(ComplexF64, Nx)
    for i in 1:Nx
        ψ_ana[i] = gaussian_wavepacket(x[i], t_last, a0, p0)
    end
    plot!(x, abs.(ψ_ana), label="Analytical", linewidth=2, dpi=150)

    savefig("IMG_orig_itime_" * string(Nt) * ".png")

end

main()