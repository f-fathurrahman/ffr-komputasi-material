using Printf

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("Atoms.jl")
include("MullerBrown.jl")
include("NEB.jl")

function plot_MullerBrown( images::Vector{Atoms}; filsave="IMG_MullerBrown_path.pdf" )
    mb = MullerBrown()
    Nx = 200
    Ny = 200
    xgrid = range(-1.7, stop=1.0, length=Nx)
    ygrid = range(-0.5, stop=2.0, length=Ny)
    Vpot = zeros(Float64,Nx,Ny)
    forces = zeros(Float64,3,1) # not used yet
    for j in 1:Ny, i in 1:Nx
        x = xgrid[i]
        y = ygrid[j]
        Vpot[i,j] = calc_energy_forces!(mb, x, y, forces)
        # Cutoff the potential (to obtain easier visualization)
        if Vpot[i,j] > 1.0
            Vpot[i,j] = 1.0
        end
    end

    plt.clf()
    plt.contour(xgrid, ygrid, transpose(Vpot), levels=10)
    plt.axis("equal")
    for i in 1:length(images)
        x = images[i].positions[1,1]
        y = images[i].positions[2,1]
        plt.plot([x], [y], marker="o", color="red")
    end
    plt.tight_layout()
    plt.savefig(filsave)
end


function create_images(;Nimages=10)
    initial = Atoms()
    initial.positions[1,1] = -0.5582247
    initial.positions[2,1] =  1.44172487

    final = Atoms()
    final.positions[1,1] = 0.62362776
    final.positions[2,1] = 0.02813212

    images = Vector{Atoms}(undef,Nimages)
    images[1] = initial
    images[Nimages] = final
    for i in 2:Nimages-1
        images[i] = Atoms()
    end
    return images
end

function main()
    images = create_images()
    neb = NEBCalculator(images)
    plot_MullerBrown(neb.images)
    println("Pass here")
end

main()