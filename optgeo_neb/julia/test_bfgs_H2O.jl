using Printf
using LinearAlgebra

include("constants.jl")
include("Atoms.jl")
include("JDFTx.jl")

function initial_hessian( Natoms )
    h = 70.0 * eV2Ha/(ANG2BOHR^2) #/(Ha2eV) * (ANG2BOHR^2) # convert to au
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
    
    NiterMax = 15
    fmax_conv = 0.01*(eV2Ha/ANG2BOHR)
    
    #println("fmax_conv = ", fmax_conv)
    #h = 70.0 * eV2Ha/(ANG2BOHR^2)
    #println("h = ", h)
    #exit()

    # we are using the same unit as ASE
#    atoms = Atoms(xyz_string="""
#    3
#
#    O   0.0  0.0   0.0
#    H   1.0  0.0   0.0
#    H   0.0  1.1   0.0
#    """,
#    LatVecs=diagm(0 => 16.0*ones(3))
#    )

    atoms = Atoms(xyz_string="""
    3

    O   0.0                   0.0                  0.0
    H   0.972786466374661    -0.1171862933028715   0.0
    H  -0.11450442319919593   0.972786466374661    0.0
    """,
    LatVecs=diagm(0 => 16.0*ones(3))
    )

    # Set up calculator
    calc = JDFTxCalculator()
    calc.Ncore = 2
    calc.use_smearing = false
    calc.kpoint_folding[:] = [1,1,1]
    calc.idx_atom_fixed = [1]
    calc.prefix_dir = "rundir_jdftx_H2O"

    Natoms = atoms.Natoms
    
    compute!(calc, atoms)

    println("Initial r  =")
    display(atoms.positions'); println()
    
    energy = atoms.energy
    forces = atoms.forces
    println("Energy = ", energy)
    println("Initial forces = ")
    display(forces'); println()
    
    MAXSTEP = 0.04*ANG2BOHR # convert to bohr

    f = vec(copy(forces))
    r = vec(copy(atoms.positions))
    r0 = zeros(size(r))
    f0 = zeros(size(f))

    H = initial_hessian(Natoms)
    
    fmax = sqrt( maximum(sum(forces.^2, dims=1)) )

    for iter = 1:NiterMax

        println("Hessian = ")
        display(H); println()

        ω, V = eigen( Symmetric(H) )
        dr = V * (V'*f ./ abs.(ω))
        steplengths = sqrt.(sum( dr.^2, dims=1 ))
        maxsteplength = maximum(steplengths)
        println("maxsteplength = ", maxsteplength)
        if maxsteplength >= MAXSTEP
            println("Scaling dr")
            dr = dr * MAXSTEP / maxsteplength
        end

        r0 = copy(r)
        f0 = copy(f)
        H_old = copy(H)

        #r[:] = r[:] + dr[:]
        # update positions
        for ia in 1:Natoms
            dr3 = reshape(dr,3,Natoms) # FIXME
            if !(ia in calc.idx_atom_fixed)
                atoms.positions[1,ia] = atoms.positions[1,ia] + dr3[1,ia]
                atoms.positions[2,ia] = atoms.positions[2,ia] + dr3[2,ia]
                atoms.positions[3,ia] = atoms.positions[3,ia] + dr3[3,ia]
            end
        end
        r = vec(copy(atoms.positions))

        energy_old = energy

        @printf("\nIter = %3d, Etot = %18.10f, fmax=%f\n", iter, energy, fmax)
        println("Forces = ")
        display(forces'); println()
        println("dr = ")
        display(reshape(dr,(3,Natoms))'); println()
        println("r  =")
        display(atoms.positions'); println()
        println("r  (in Angstrom) =")
        display(atoms.positions'*BOHR2ANG); println()

        compute!(calc, atoms)
        energy = atoms.energy
        forces[:] = atoms.forces[:]
        
        fmax = sqrt( maximum(sum(forces.^2, dims=1)) )
        if fmax < fmax_conv
            println("BFGS is converged")
            break
        end

        f = vec(copy(forces))
        H = update_hessian( H_old, r, f, r0, f0 )

    end


end

main()