using Dierckx: Spline1D, derivative
using Printf


function initial_conditions( Np, x0, β, m, E0; xmin=-0.6, xmax=0.6)

    # Particle spacings
    h = (xmax - xmin)/(Np-1)

    # Normalization for an initial Gaussian wave function
    anorm = (2.0*β/π)^(1.0/4.0)
    
    #Particle positions
    x = zeros(Np)
    
    # Build initial wave packet amplitude and density
    R = zeros(Np)
    C = zeros(Np)
    ρ = zeros(Np)
    for i in 1:Np
        x[i] = xmin + (i-1)*h
        R[i] = anorm*exp( -β*( ( x[i] - x0)^2 ) )
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
    x0 = 0.0
    
    # width parameter for Gaussian wavepacket
    β = 9.0
    
    # Initial translational energy
    E0 = 0.0
        
    # Mass in a.u.
    m = 1.0  # electron

    Np = 100

    x, R, C, ρ, vx, delv, phase = initial_conditions( Np, x0, β, m, E0, xmin=-0.1, xmax=0.1 )


    d1x = zeros(Np)
    d2x = zeros(Np)
    pot = zeros(Np)
    Q = zeros(Np)
    pot = zeros(Np)

    Ntime = 10
    
    f = open("TEMP_trajectory.dat", "w")

    for k = 1:Ntime
    
        t = k*dt

        println("t = ", t)
    
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
    
        #plt.plot( t*ones(Np), x, marker="o", linewidth=0, color="black" )

        @printf(f, "%18.10f ", t)
        for i in 1:Np
            @printf(f, "%18.10f ", x[i])
        end
        @printf(f, "\n")
    
    end

    #plt.savefig("TEMP_traj_x.png", dpi=150)

    close(f)

end

main()