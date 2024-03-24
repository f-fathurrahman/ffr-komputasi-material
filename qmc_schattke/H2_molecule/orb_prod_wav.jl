function orb_prod_wav!(
    global_vars, psi
)

    Norbitals = size(psi, 1)
    @assert Norbitals == 1

    R_electrons = global_vars.R_electrons
    R_atoms = global_vars.R_atoms
    wavec = global_vars.wavec
    α = global_vars.α
    Natoms = global_vars.Natoms

    phi = zeros(Float64, Norbitals, Natoms)
    dr = zeros(Float64, 3)

    for ia in 1:Natoms
        dr[1:3] .= R_electrons[1:3] .- R_atoms[1:3,ia]
        s = sqrt( sum(dr[1:3].^2) ) 
        phi[1,ia] = (1.0 + wavec*s)*exp(-α*s)
    end
    # Associate single particle wavefunction with electron, here 1 to 1
    psi[1,1] = phi[1,1]*phi[1,2]
    psi[1,2] = phi[1,1]*phi[1,2]
    return
end




