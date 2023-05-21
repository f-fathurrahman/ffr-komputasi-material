using LinearAlgebra
using Printf

using Plots
using PlotThemes
theme(:dark)

const METER_TO_BOHR = 18897259885.789

# [1D] integration - Simpson's 1/3 rule
#      f function    a = lower bound    b = upper bound
#      Must have odd number of data points
#      Simpson's coefficients   1 4 2 4 ... 2 4 1

function simpson1d(f, a, b)
    numS = length(f) # number of data points

    @assert mod(numS,2) == 1

    sc = 2*ones(numS)
    sc[2:2:numS-1] .= 4
    sc[1] = 1
    sc[numS] = 1

    h = (b - a)/(numS - 1)

    return (h/3) * sum(f .* sc)
end


Nx = 1001 # must be an odd number - number of grid points
Nt = 10000 # number of time steps

me = 1.0  # electron mass
hbar = 1.0 # hbar Planck's constant
e = 1.0  # elementary charge

L = 4e-9 * METER_TO_BOHR

Vpot = zeros(Float64, Nx)

ψR = zeros(Float64, Nx)
ψI = zeros(Float64, Nx)

x0 = L/4 # pulse center
s = L/25
wL = 1.6e-10 * METER_TO_BOHR

h = hbar*2π
KE = (h/wL)^2/(2*me*e) # theoretical KE
x = range(0, stop=L, length=Nx)
dx = x[2] - x[1];
t = zeros(Float64, Nt)
C1 = 1/10 # stability factor
dt = C1 * 2 * me * dx^2 / hbar
C2 = e*dt/hbar
C3 = -hbar^2 / (2 * me * dx^2 * e)
T = Nt * dt

# Initial wave packet
for i in 1:Nx
    gauss_f = exp(-0.5*((x[i] - x0)/s)^2)
    ψR[i] = gauss_f * cos( 2π*(x[i] - x0)/wL )
    ψI[i] = gauss_f * sin( 2π*(x[i] - x0)/wL )
end

# Normalize initial wave packet
ψ2 = ψR.^2 + ψI.^2
A = simpson1d(ψ2, 0.0, L)
#A = sum(ψ2)*dx
ψR = ψR / sqrt(A)
ψI = ψI / sqrt(A)

str_itime = @sprintf("%06d", 0)
plot(x, ψR, label="Re", linewidth=2, dpi=150)
ylims!(-0.5, 0.5)
savefig("IMG_psi_R_" * str_itime * ".png")

plot(x, ψI, label="Im", linewidth=2, dpi=150)
ylims!(-0.5, 0.5)
savefig("IMG_psi_I_" * str_itime * ".png")

prob_density = ψR.^2 + ψI.^2
@printf("Check norm: %18.10f\n", sum(prob_density)*dx)
ψR1 = copy(ψR)
ψI1 = copy(ψI)
ψP11 = copy(prob_density)

# Kinetic energy KE
fn = zeros(ComplexF64, Nx-1)
for i in 2:Nx-1
    fn[i] = C3 * ( ψR[i] - im*ψI[i] ) * 
            (ψR[i+1] - 2*ψR[i] + ψR[i-1] + im*( ψI[i+1] - 2*ψI[i] + ψI[i-1]) )
end
K1avg = simpson1d(fn[2:end], 0, L)
println("K1avg = ", K1avg)

#
# FDTD solution of Schrodinger Equation
#
cs = 1.0
probD = zeros(Float64, Nx, 9)

for itime in 1:Nt
    for i in 2:Nx-1
        ψR[i] = ψR[i] - C1*( ψI[i+1] - 2*ψI[i] + ψI[i-1] ) + C2*Vpot[i]*ψI[i]
    end
    for i in 2:Nx-1
        ψI[i] = ψI[i] + C1*( ψR[i+1] - 2*ψR[i] + ψR[i-1] ) - C2*Vpot[i]*ψR[i]
    end
    if itime % 50 == 0
        #
        println("itime = ", itime, ", do plot")
        #
        plot(x, ψR, label="Re", linewidth=2, dpi=150)
        ylims!(-0.5, 0.5)
        str_itime = @sprintf("%06d", itime)
        savefig("IMG_psi_R_" * str_itime * ".png")
        #
        plot(x, ψI, label="Im", linewidth=2, dpi=150)
        ylims!(-0.5, 0.5)
        savefig("IMG_psi_I_" * str_itime * ".png")
    end
end