using Printf

import PyPlot
const plt = PyPlot

plt.rc("text", usetex=true)

include("Atoms.jl")
include("MullerBrown.jl")

function do_plot()
    mb = MullerBrown()
    Nx = 200
    Ny = 200
    xgrid = range(-1.7, stop=0.5, length=Nx)
    ygrid = range(-0.5, stop=2.0, length=Ny)
    Vpot = zeros(Float64,Nx,Ny)
    forces = zeros(Float64,3,1) # not used
    for j in 1:Ny, i in 1:Nx
        x = xgrid[i]
        y = ygrid[j]
        Vpot[i,j] = calc_energy_forces!(mb, x, y, forces)
    end

    plt.clf()
    #plt.surf(xgrid, ygrid, transpose(Vpot))
    #plt.contourf(xgrid, ygrid, transpose(Vpot), levels=20)
    plt.contour(xgrid, ygrid, transpose(Vpot), levels=20)
    plt.axis("square")
    plt.axis("equal")
    #plt.xticks([]); plt.yticks([])
    #plt.text(-4.5, 4.5, "State "*string(ist))
    plt.tight_layout()
    plt.savefig("IMG_MullerBrown.pdf")
end

do_plot()