using Dierckx: Spline1D, derivative
using Printf
import PyPlot
const plt = PyPlot

function initial_conditions( Np, x1, x2, β, m, E0; xmin=-0.6, xmax=0.6)

    @assert Np % 2 == 0

    # Particle spacings
    h = (xmax - xmin)/(Np-1)
    println("h = ", h)

    # Normalization for an initial Gaussian wave function
    anorm = (2.0*β/π)^(1.0/4.0)
    
    #Particle positions
    x = zeros(Np)

    Np2 = round(Int64,Np/2)
    xtmp = zeros(Np2)
    for i in 1:Np2
        xtmp[i] = -0.01 + (i-1)*0.02/(Np2-1)
    end
    x[1:Np2] = -0.1 .+ xtmp[:]
    x[Np2+1:end] = 0.1 .+ xtmp[:]
    
    # Build initial wave packet amplitude and density
    R = zeros(Np)
    C = zeros(Np)
    ρ = zeros(Np)
    for i in 1:Np
        R[i] = 0.5 * anorm * ( exp(-β*((x[i] - x1)^2)) + exp(-β*((x[i] - x2)^2)) )
        #R[i] = 0.5 * anorm * exp(-β*( (x[i] - x1)^2 ))
        C[i] = log(R[i])
        ρ[i] = R[i]^2
    end
    
    vx = zeros(Np)
    delv = zeros(Np)
    # Initial particle velocities (all particles are initially the 
    # same velocity and 0 divergence)
    # These velocities are obtained from the equation
    # Etrans = 1/2 * m * v^2 .
    for i in 1:Np
        vx[i] = sqrt(2.0*E0/m)
        delv[i] = 0.0
    end
    
    phase = zeros(Np)
    # Initial phase (action function) for each particle
    for i in 1:Np
        phase[i] = sqrt(2.0*m*E0)*x[i]
    end
    
    return x, R, C, ρ, vx, delv, phase
end


function calc_deriv!( x, fx, d1x, d2x )
    # interpolant for fx
    spl_fx = Spline1D(x, fx)
    d1x[:] = derivative(spl_fx, x)
    spl_d1x = Spline1D(x, d1x)
    d2x[:] = derivative(spl_d1x, x)
    return
end


function main()

    # timestep in a.u.
    dt = 1e-2
    
    # center of Gaussian function
    x1 = -0.1

    x2 = 0.1
    
    # width parameter for Gaussian wavepacket
    β = 1000.0
    
    # Initial translational energy
    E0 = 0.0
        
    # Mass in a.u.
    m = 1.0  # electron

    Np = 10

    x, R, C, ρ, vx, delv, phase = initial_conditions( Np, x1, x2, β, m, E0, xmin=-0.25, xmax=0.25 )

    #plt.clf()
    #plt.plot(x, R)
    #plt.savefig("TEMP_initial.png")


    d1x = zeros(Np)
    d2x = zeros(Np)
    pot = zeros(Np)
    Q = zeros(Np)
    pot = zeros(Np)

    Ntime = 5
    
    f = open("TEMP_trajectory.dat", "w")

    for k = 1:Ntime
    
        t = k*dt

        @printf("t = %18.10f\n", t)
    
        calc_deriv!(x, C, d1x, d2x)
    
        # Calculate quantum and classical potentials
        for i in 1:Np
            # Quantum Potential.
            Q[i] = -1.0/(2.0*m)*(d2x[i] + d1x[i]^2)
            # Classical Eckart Potential.
            #pot[i] = vb*sech(wx*(x[i] - x0))^2
            # Harmonic potential
            #pot[i] = 0.5*x[i]^2
        end
    
        # Calculate phase (S) using potentials (Quantum Lagrangian: T - (V + Q))
        for i in 1:Np
            KE = 0.5*m*vx[i]^2
            lagrange = KE - ( pot[i] + Q[i] )
            phase[i] = phase[i] + lagrange*dt
        end
    
        # Update particle positions
        for i in 1:Np
            x[i] = x[i] + vx[i]*dt
        end
    
        idx_sorted = sortperm(x)
        x[:]     = x[idx_sorted]
        phase[:] = phase[idx_sorted]
        #vx[:]    = vx[idx_sorted]
        #delv[:]  = delv[idx_sorted]
        ρ[:]     = ρ[idx_sorted]
        C[:]     = C[idx_sorted]
    
        # Call subroutine for spatial derivatives of phase
        calc_deriv!(x, phase, d1x, d2x)
    
        # Update velocities and probability density
        for i in 1:Np
            vx[i] = (1.0/m)*d1x[i]
            delv[i] = (1.0/m)*d2x[i]
            ρ[i] = ρ[i]*exp( -delv[i]*dt )
        end
    
        @printf(f, "%18.10f ", t)
        for i in 1:Np
            @printf(f, "%18.10f ", x[i])
        end
        @printf(f, "\n")
    
    end

    close(f)


end

main()