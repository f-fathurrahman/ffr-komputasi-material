
mutable struct GlobalVariables
    Natoms::Int64
    Nelectrons::Int64
    NelectronsPerSpin::Int64
    Norbitals::Int64
    R_atoms::Matrix{Float64}
    R_electrons::Matrix{Float64}
    R_electrons_new::Matrix{Float64}
    DOLD::Vector{Float64}
    DNEW::Vector{Float64}
    β1::Float64
    β2::Float64
    γ::Float64
    CJAS::Float64
    α::Float64
    wavec::Float64
end

include("orb_prod_wav.jl")

function create_global_vars()

    Natoms = 2
    Nelectrons = 2
    NelectronsPerSpin = 1
    Norbitals = 1

    distance_H2 = 1.40 # distance of both H2 nuclei
    R_atoms = zeros(Float64, 3, Natoms)
    R_atoms[1,2] = distance_H2 # x coord of 2nd atom

    R_electrons = zeros(Float64, 3, Nelectrons)
    R_electrons_new = zeros(Float64, 3, Nelectrons)

    # For Jastrow factor evalution ???
    DOLD = zeros(Float64, Nelectrons)
    DNEW = zeros(Float64, Nelectrons)

    # Jastrow
    β1 = 0.01
    β2 = 0.02
    γ = 0.0001
    CJAS = 0.00001

    α = 0.679
    wavec = 0.1

    return GlobalVariables(
        Natoms, Nelectrons, NelectronsPerSpin, Norbitals,
        R_atoms, R_electrons, R_electrons_new,
        DOLD, DNEW,
        β1, β2,
        γ, CJAS, α, wavec
    )
end


function initial_electron_positions!(global_vars)
    
    Nelectrons = global_vars.Nelectrons
    R_electrons = global_vars.R_electrons
    R_electrons_new = global_vars.R_electrons_new
    Norbitals = global_vars.Norbitals
    NelectronsPerSpin = global_vars.NelectronsPerSpin
    DOLD = global_vars.DOLD

    Hpsi = zeros(Float64, Norbitals, Nelectrons)

    iorbital = 1 # only one orbital

    # Random initial electron positions
    for iel in 1:Nelectrons
        ispin = 1 # what's this?
        if iel > NelectronsPerSpin
            ispin = 2
        end
        for i in 1:3
            rd = rand(Float64) - 0.5
            R_electrons[i,iel] = R_electrons[i,iel] + rd   # relative from nucleus
            R_electrons_new[i,iel] = R_electrons[i,iel]  # also assign RNEU <= RE
        end
        orb_prod_wav!( global_vars, Hpsi ) # evaluate psi at this position
        DOLD[ispin] = Hpsi[iorbital,iel]
    end

    return
end


#function debug_main()


    global_vars = create_global_vars()

    LENGTH = 10.0 # size of arb. box to display positions

    # For debugging purpose
    MCPRE = 1000   # 100000
    MCMAX = 20000  # 2000000

    NDIV = 21 # NDIV must be odd

    # Maximum step width
    STEPMAX = 1.0

    initial_electron_positions!(global_vars)

#    return
#end
#debug_main()
