# 1d FDTD simulation of psi
# Sullivan, Quantum Mechanics for Electrical Engineers

using Printf
using LinearAlgebra

import PyPlot
const plt = PyPlot

# FIXME: Use atomic units?

function main()

# Number of points in the problem space.
NN = 400

# Planck constant
ħ = 1.054e-34

# Free space mass of an electron
m0 = 9.1e-31

# Effective mass: Si is 1.08, Ge is 0.067, GaAs is 0.55
meff = 1.0

# Mass of an electron
melec = meff*m0

# Charge of an electron
ecoul = 1.6e-19

# Dielectric of free space
epsz = 8.85e-9

# Energy conversion factors
eV2J = 1.6e-19
J2eV = 1/eV2J

# The cell size
Δx = 0.1e-9 # m

# Time steps
Δt = 2e-17 # s

# Ra must be < 0.15
Ra = (0.5*ħ/melec)*(Δt/Δx^2)

println("Ra = ", Ra)

# Specify the potential
V = zeros(Float64, NN)

NN_half = Int64(NN/2)
# Barrier
#for n  in NN_half:(NN_half+50)
#    V[n] = 0.15*eV2J
#end

# Semiconductor conduction band
for n in 1:NN_half
    V[n] = 0.1*eV2J
end

#for n=NN/2+1:NN
#    %V(n) = .2*eV2J
#end

# Electric ﬁeld
#for n=1:NN
#    V(n) = -(0.2/400)*(n-400)*eV2J
#end

# Initialize a sine wave in a gaussian envelope

# Pulse wavelength
λ = 50
# lambda = 25

σ = 50 # Pulse width

nc = 150 # Starting position

prl = zeros(Float64, NN) # The real part of the state variable
pim = zeros(Float64, NN) # The imaginary part of the state variable
ptot = 0.0
for n in 2:NN-1
    prl[n] = exp(-1.0*((n-nc)/σ)^2)*cos(2*pi*(n-nc)/λ) 
    pim[n] = exp(-1.0*((n-nc)/σ)^2)*sin(2*pi*(n-nc)/λ) 
    ptot = ptot + prl[n]^2 + pim[n]^2
end
pnorm = sqrt(ptot) # Normalization constant

# Normalize and check
ptot = 0.0
for n in 1:NN
    prl[n] = prl[n]/pnorm
    pim[n] = pim[n]/pnorm
    ptot = ptot + prl[n]^2 + pim[n]^2
end

@printf("ptot = %f (should be 1.0)\n", ptot) # This should have the value 1


T = 0
n_step = 1000


# Cell size in nm.
dx_nm = Δx*1e9
println("Cell size (nm) = ", dx_nm)

# Length in nm for plotting
xgrid = dx_nm:dx_nm:dx_nm*NN

psi2 = zeros(Float64, NN)

# This is the core FDTD program
for m in 1:n_step
    T = T + 1
    for n in 2:NN-1
        prl[n] = prl[n] - Ra*( pim[n-1] -2*pim[n] + pim[n+1] ) + (Δt/ħ)*V[n]*pim[n]
    end
    for n in 2:NN-1
        pim[n] = pim[n] + Ra*( prl[n-1] -2*prl[n] + prl[n+1] ) - (Δt/ħ)*V[n]*prl[n]
    end
    for n in 1:NN
        psi2[n] = prl[n]^2 + pim[n]^2
    end
    plt.clf()
    plt.plot(xgrid, psi2, label="psi2")
    plt.savefig("IMG_Se1_1_psi_"*string(m)*".png", dpi=150)
    @printf("Step %d is done\n", m)
end

# ------------------------------------------------
# Calculate the expected values
PE = 0.0
psi = zeros(ComplexF64, NN)
for n in 1:NN
    # Write as a complex function
    psi[n] = prl[n] + im*pim[n]
    PE = PE + psi[n]*conj(psi[n])*V[n]
end
    
# This checks normalization
println("Check norm = ", dot(psi, psi))
    
PE = PE*J2eV # Potential energy

ke = 0.0 + im*0.0
for n in 2:NN-1
    lap_p = psi[n+1] - 2*psi[n] + psi[n-1]
    ke = ke + lap_p*conj(psi[n])
end
# Kinetic energy
KE = -J2eV*( (ħ/Δx)^2 / (2*melec) )*real(ke)

end

main()
