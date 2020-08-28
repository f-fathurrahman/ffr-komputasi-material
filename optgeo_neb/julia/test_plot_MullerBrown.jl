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
    xgrid = range(-1.7, stop=1.0, length=Nx)
    ygrid = range(-0.5, stop=2.0, length=Ny)
    Vpot = zeros(Float64,Nx,Ny)
    forces = zeros(Float64,3,1) # not used
    for j in 1:Ny, i in 1:Nx
        x = xgrid[i]
        y = ygrid[j]
        Vpot[i,j] = calc_energy_forces!(mb, x, y, forces)
        # Cutoff the potential (to obtain easier visualization)
        if Vpot[i,j] > 1.0
            Vpot[i,j] = 1.0
        end
    end

    println("maximum(Vpot) = ", maximum(Vpot))

    x1 = -0.5582247
    y1 =  1.44172487

    x2 = 0.62362776
    y2 = 0.02813212

    plt.clf()
    #plt.surf(xgrid, ygrid, transpose(Vpot))
    #plt.contourf(xgrid, ygrid, transpose(Vpot), levels=20)
    plt.contour(xgrid, ygrid, transpose(Vpot), levels=10)
    plt.axis("equal")
    plt.plot([x1], [y1], marker="*", color="black")
    plt.plot([x2], [y2], marker="o", color="red")
    plt.tight_layout()
    plt.savefig("IMG_MullerBrown.pdf")
end

do_plot()