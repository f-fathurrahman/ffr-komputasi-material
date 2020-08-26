using Printf
using LinearAlgebra


include("Atoms.jl")
include("MullerBrown.jl")

function initial_hessian( Natoms )
    h = 70.0 #/(2*Ry2eV) * (ANG2BOHR^2)
    #h = 70.0
    H = diagm( 0 => h*ones(3*Natoms) )
    return H
end
    
function update_hessian( H_old, r, f, r0, f0 )
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
    
    # we are using the same unit as ASE
    atoms = Atoms(1)
    atoms.positions[1,1] = -0.55
    atoms.positions[2,1] =  1.30

    mb = MullerBrown()

    Natoms = atoms.Natoms
    
    calc_energy_forces!(mb, atoms)

    println("Initial r  =")
    display(atoms.positions'); println()
    
    energy = atoms.energy
    forces = atoms.forces
    println("Energy = ", energy)
    println("Initial forces = ")
    display(forces'); println()
    
    MAXSTEP = 0.04 # !! in angstrom

    f = vec(copy(forces))
    r = vec(copy(atoms.positions))
    r0 = zeros(size(r))
    f0 = zeros(size(f))

    H = initial_hessian(Natoms)

    NiterMax = 15
    for iter = 1:NiterMax

        println("Hessian = ")
        display(H); println()

        omega, V = eigen( Symmetric(H) )
        dr = V * (V'*f ./ abs.(omega))
        steplengths = sqrt.(sum( dr.^2, dims=1 ))
        maxsteplength = maximum(steplengths)
        if maxsteplength >= MAXSTEP
            println("Scaling dr")
            dr = dr * MAXSTEP / maxsteplength
        end

        r0 = copy(r)
        f0 = copy(f)
        H_old = copy(H)

        r[:] = r[:] + dr[:]
        
        # update positions
        r_new = reshape(r,3,Natoms)
        atoms.positions[1,1] = r_new[1,1]
        atoms.positions[2,1] = r_new[2,1]

        energy_old = energy

        @printf("\nIter = %3d, Etot = %18.10f\n", iter, energy)
        println("Forces = ")
        display(forces'); println()
        println("dr = ")
        display(reshape(dr,(3,Natoms))'); println()
        println("r  =")
        display(atoms.positions'); println()

        calc_energy_forces!(mb, atoms)
        energy = atoms.energy
        forces[:] = atoms.forces[:]
        
        fmax = norm(forces)
        if fmax < 0.01
            println("BFGS is converged")
            break
        end

        f = vec(copy(forces))
        H = update_hessian( H_old, r, f, r0, f0 )

    end


end

main()