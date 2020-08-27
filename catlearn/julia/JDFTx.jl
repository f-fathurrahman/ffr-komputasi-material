# Contain options for JDFTx calculation
mutable struct JDFTxCalculator
    exe_command::String
    input_file::String
    output_file::String
    prefix_dir::String
    clean_dir::Bool
    kpoint_folding::Array{Int64}
    idx_atom_fixed::Array{Int64}
    Ncore::Int64
    use_smearing::Bool
end

function JDFTxCalculator()
    exe_command = "jdftx-1.6.0 -c 1 -i"
    prefix_dir = "rundir_jdftx"
    input_file = "inp"
    output_file = "out"
    clean_dir = false
    kpoint_folding = [1,1,1]
    idx_atom_fixed = [0]
    Ncore = 1
    use_smearing = false
    return JDFTxCalculator(exe_command, input_file, output_file, prefix_dir,
        clean_dir, kpoint_folding, idx_atom_fixed, Ncore, use_smearing)
end

function write_jdftx(calc::JDFTxCalculator, atoms::Atoms)

    input_file = calc.input_file
    prefix_dir = calc.prefix_dir
    use_smearing = calc.use_smearing

    if !Base.Filesystem.isdir(prefix_dir)
        Base.Filesystem.mkdir(prefix_dir)
    end
    f = open(prefix_dir*"/$input_file", "w")
    
    # Unit cell
    LatVecs = atoms.LatVecs
    @printf(f, "lattice \\\n")
    @printf(f, "%18.10f %18.10f %18.10f \\\n", LatVecs[1,1], LatVecs[1,2], LatVecs[1,3])
    @printf(f, "%18.10f %18.10f %18.10f \\\n", LatVecs[2,1], LatVecs[2,2], LatVecs[2,3])
    @printf(f, "%18.10f %18.10f %18.10f \n",   LatVecs[3,1], LatVecs[3,2], LatVecs[3,3])
    @printf(f, "\n")

    # Species and pseudopotentials
    @printf(f, "ion-species GBRV/\$ID_pbe.uspp\n")

    # Cutoff (default)
    @printf(f, "elec-cutoff 20 100\n")

    # Ionic positions
    atpos = atoms.positions
    atsymbs = atoms.atsymbs
    idx_atom_fixed = calc.idx_atom_fixed
    @printf(f, "coords-type cartesian\n")
    for ia in 1:atoms.Natoms
        if ia in idx_atom_fixed
            @printf(f, "ion %s %18.10f %18.10f %18.10f 0\n",
                atsymbs[ia], atpos[1,ia], atpos[2,ia], atpos[3,ia])
        else
            @printf(f, "ion %s %18.10f %18.10f %18.10f 1\n",
                atsymbs[ia], atpos[1,ia], atpos[2,ia], atpos[3,ia])
        end
    end

    kfold = calc.kpoint_folding
    @printf(f, "kpoint-folding %d %d %d\n", kfold[1], kfold[2], kfold[3])
    
    if use_smearing
        # Default
        @printf(f, "elec-smearing cold 0.01\n")
    end

    @printf(f, "dump-name \$VAR\n")
    @printf(f, "initial-state \$VAR\n")

    # For isolated molecule
    #coulomb-interaction isolated
    #coulomb-truncation-embed 8 8 8.00001
    
    @printf(f, "dump End State\n")
    @printf(f, "dump End Forces\n")
    @printf(f, "dump End Ecomponents\n")

    close(f)
end

function read_energy(calc::JDFTxCalculator)
    f = open(joinpath(calc.prefix_dir, "Ecomponents"))
    energy = 0.0
    while !eof(f)
        line = readline(f)
        if occursin("Etot =", line)
            ll = split(line, "=", keepempty=false)
            energy = parse(Float64,ll[2])
        end
        # Prefer F rather that Etot in case of use_smearing=true
        if occursin("F =", line)
            ll = split(line, "=", keepempty=false)
            energy = parse(Float64,ll[2])
        end
    end
    return energy
end


function read_forces!(calc::JDFTxCalculator, forces::Array{Float64})
    Natoms = size(forces,2)
    f = open(joinpath(calc.prefix_dir, "force"))
    ia = 0
    while !eof(f)
        line = readline(f)
        if occursin("#", line) continue end
        if !occursin("force", line) continue end
        ll = split(line, keepempty=false)
        ia = ia + 1
        if ia > Natoms
            error("Error reading force, too much data")
        end
        forces[1,ia] = parse(Float64, ll[3])
        forces[2,ia] = parse(Float64, ll[4])
        forces[3,ia] = parse(Float64, ll[5])
    end
    close(f)
end

# return energy, overwrite forces
function compute!(
    calc::JDFTxCalculator, atoms::Atoms, forces::Array{Float64,2}
)

    prefix_dir = calc.prefix_dir
    input_file = calc.input_file

    if calc.clean_dir
        run(`rm -rfv $prefix_dir`)
    end
    
    write_jdftx(calc, atoms)

    Ncore = calc.Ncore
    output_file = calc.output_file
    # assume relative path
    cd("./$prefix_dir")
    run(pipeline(`jdftx-1.6.0 -c $Ncore -i $input_file`, stdout="$output_file"))
    cd("../")
    # need to change this?, use pwd to save the original directory
    # not needed right now as we assume relative path
    
    energy = read_energy(calc)
    read_forces!(calc, forces)

    return energy
end

