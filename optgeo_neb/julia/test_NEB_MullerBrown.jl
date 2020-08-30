using Printf
using LinearAlgebra

import PyPlot
const plt = PyPlot
plt.rc("text", usetex=true)

include("Atoms.jl")
include("MullerBrown.jl")
include("NEB.jl")

function plot_MullerBrown(
    mb::MullerBrown, images::Vector{Atoms};
    filsave="IMG_MullerBrown_path.pdf"
)
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

function create_images(;Nimages=12)
    initial = Atoms()
    initial.positions[1,1] = -0.5582247033
    initial.positions[2,1] =  1.4417248736

    final = Atoms()
    final.positions[1,1] = 0.6236277588
    final.positions[2,1] = 0.0281321159

    images = Vector{Atoms}(undef,Nimages)
    images[1] = initial
    images[Nimages] = final
    for i in 2:Nimages-1
        images[i] = Atoms()
    end
    return images
end

function initial_hessian(Ndofs)
    h = 70.0
    H = diagm( 0 => h*ones(Ndofs) )
    return H
end
    
function update_hessian(H_old, r, f, r0, f0)
    Natoms = size(r,2)
    dr = r - r0
    # FIXME: need this?
    if maximum( abs.(dr) ) < 1e-7
        return diagm( 0 => 70*ones(3*Natoms) )
    end
    df = f - f0
    a = dot(dr, df)
    dg = H_old*dr
    b = dot(dr, dg)
    return H_old - (df*df')/a - (dg*dg')/b
end

function main()
    images = create_images(Nimages=10)
    neb = NEBCalculator(images)
    Nimages = neb.Nimages
    Natoms = 1

    mb = MullerBrown()
    #plot_MullerBrown(mb, neb.images)

    setup_initial_final!(mb, neb)
    println("ene initial = ", neb.images[1].energy)
    println("ene final   = ", neb.images[end].energy)
    compute!(mb, neb)
    
    println("Initial energies:")
    for i in 1:Nimages
        x = neb.images[i].positions[1,1]
        y = neb.images[i].positions[2,1]
        @printf("%3d r=[%18.10f,%18.10f] E=%18.10f\n", i, x, y, images[i].energy)
    end

    forces = neb.neb_forces
    
    MAXSTEP = 0.04 # !! in angstrom

    # using linear indexing
    f = vec(copy(forces))
    r = vec(copy(get_moving_positions(neb)))

    H = initial_hessian(3*Natoms*(Nimages-2))

    # Calculate fmax (FIXME: check against ASE implementation)
    ff = reshape(forces, (3,Natoms*(Nimages-2)))
    fmax = sqrt( maximum(sum(ff.^2, dims=1)) )

    NiterMax = 40
    for iter = 1:NiterMax

        #println("Hessian = ")
        #display(H); println()

        omega, V = eigen( Symmetric(H) )
        println("omega = ", omega)
        dr = V * (V'*f ./ abs.(omega))
        dr = reshape(dr, 3,Natoms*(Nimages-2)) # FIXME: remove this?
        println("dr = ", dr)
        steplengths = sqrt.(sum( dr.^2, dims=1 ))
        println("MyBFGS: steplengths = ", steplengths)
        maxsteplength = maximum(steplengths)
        if maxsteplength >= MAXSTEP
            println("Scaling dr")
            println("maxsteplength = ", maxsteplength)
            dr = dr * MAXSTEP / maxsteplength
        end

        # Before update positions,copy old values
        r0 = copy(r)
        f0 = copy(f)
        H_old = copy(H)

        r[:] = r[:] + dr[:]
        
        # update positions
        r_new = reshape(r, (3,Natoms,Nimages-2))
        for i in 2:Nimages-1
            images[i].positions[1,1] = r_new[1,1,i-1]
            images[i].positions[2,1] = r_new[2,1,i-1]
        end

        println("New positions after update positions")
        for i in 1:Nimages
            x = neb.images[i].positions[1,1]
            y = neb.images[i].positions[2,1]
            @printf("%3d r=[%18.10f,%18.10f]\n", i, x, y)
        end

        #energy_old = energy

        #@printf("\nIter = %3d, Etot = %18.10f, fmax=%18.10f\n", iter, energy, fmax)
        @printf("\nIter = %3d, fmax=%18.10f\n", iter, fmax)

        compute!(mb, neb)
        forces[:] = neb.neb_forces[:]
        # Evaluate convergence
        ff = reshape(forces, (3,Natoms*(Nimages-2)))
        fmax = sqrt( maximum(sum(ff.^2, dims=1)) )
        if fmax < 0.05
            println("BFGS is converged")
            @printf("\nIter = %3d, fmax=%18.10f\n", iter+1, fmax)
            break
        end
        f = vec(copy(forces)) # linear index
        H = update_hessian( H_old, r, f, r0, f0 )
    end

    for i in 1:Nimages
        @printf("%3d %18.10f\n", i, images[i].energy)
    end
    plot_MullerBrown(mb, neb.images, filsave="IMG_MullerBrown_path_v2.pdf")

end

main()