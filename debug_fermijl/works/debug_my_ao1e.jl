# only one BasisSet
function my_ao1e(BS::BasisSet, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_kin_sph!
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_nuc_sph!
    else
        throw(ArgumentError("Invalid one-eletron integral name: $compute"))
    end

    # Pre allocate output
    out = zeros(T, BS.nbas, BS.nbas)

    # Pre compute a list of angular momentum numbers (l) for each shell
    lvals = [Libcint.CINTcgtos_spheric(i-1, BS.lc_bas) for i = 1:BS.nshells]
    Lmax = maximum(lvals)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset = [sum(lvals[1:(i-1)]) for i = 1:BS.nshells]

    buf = zeros(Cdouble, Lmax^2)
    # ffr: parallelizations are removed

    for i in 1:BS.nshells
        Li = lvals[i]
        ioff = ao_offset[i]
        for j in i:BS.nshells
            Lj = lvals[j]
            joff = ao_offset[j]
            # Call libcint
            libcint_1e!(buf, [i,j], BS)
            # Loop through shell block and save unique elements
            for js = 1:Lj
                J = joff + js
                for is = 1:Li
                    I = ioff + is
                    J < I ? break : nothing
                    out[I,J] = buf[is + Li*(js-1)]
                    out[J,I] = out[I,J]
                end
            end
        end
    end
    return out
end


# Two different basis set
function my_ao1e(BS1::BasisSet, BS2::BasisSet, compute::String, T::DataType = Float64)

    if compute == "overlap"
        libcint_1e! =  cint1e_ovlp_sph!
    elseif compute == "kinetic"
        libcint_1e! =  cint1e_kin_sph!
    elseif compute == "nuclear"
        libcint_1e! =  cint1e_nuc_sph!
    else
        throw(ArgumentError("Invalid one-eletron integral name: $compute"))
    end

    ATM_SLOTS = 6
    BAS_SLOTS = 8

    natm = BS1.natoms + BS2.natoms
    nbas = BS1.nbas + BS2.nbas
    nshells = BS1.nshells + BS2.nshells

    lc_atm = zeros(Cint, natm*ATM_SLOTS)
    lc_bas = zeros(Cint, nshells*BAS_SLOTS)
    env = zeros(Cdouble, length(BS1.lc_env) + length(BS2.lc_env) -20)

    # Prepare the lc_atom input 
    off = 20
    ib = 0 
    atoms = vcat(BS1.atoms, BS2.atoms)
    for i = eachindex(atoms)
        A = atoms[i]
        # lc_atom has ATM_SLOTS (6) "spaces" for each atom
        # The first one (Z_INDEX) is the atomic number
        lc_atm[1 + ATM_SLOTS*(i-1)] = Cint(A.Z)
        # The second one is the env index address for xyz
        lc_atm[2 + ATM_SLOTS*(i-1)] = off
        env[off+1:off+3] .= A.xyz ./ MyMolecules.bohr_to_angstrom
        off += 4 # Skip an extra slot reserved for nuclear model
        # The remaining 4 slots are zero.
    end

    for i = eachindex(BS1.atoms)
        # Prepare the lc_bas input
        for j = eachindex(BS1.basis[i])
            B = BS1[i,j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end

    for i = eachindex(BS2.atoms)
        # Prepare the lc_bas input
        for j = eachindex(BS2.basis[i])
            B = BS2[i,j] 
            Ne = length(B.exp)
            Nc = length(B.coef)
            # lc_bas has BAS_SLOTS for each basis set
            # The first one is the index of the atom starting from 0
            lc_bas[1 + BAS_SLOTS*ib] = i-1
            # The second one is the angular momentum
            lc_bas[2 + BAS_SLOTS*ib] = B.l
            # The third is the number of primitive functions
            lc_bas[3 + BAS_SLOTS*ib] = Nc
            # The fourth is the number of contracted functions
            lc_bas[4 + BAS_SLOTS*ib] = 1
            # The fifth is a κ parameter
            lc_bas[5 + BAS_SLOTS*ib] = 0
            # Sixth is the env index address for exponents
            lc_bas[6 + BAS_SLOTS*ib] = off
            env[off+1:off+Ne] .= B.exp
            off += Ne
            # Seventh is the env index address for contraction coeff
            lc_bas[7 + BAS_SLOTS*ib] = off
            env[off+1:off+Nc] .= B.coef
            off += Nc
            # Eigth, nothing
            ib += 1
        end
    end

    # Allocate output array
    out = zeros(T, BS1.nbas, BS2.nbas)
    lvals1 = [Libcint.CINTcgtos_spheric(i-1, BS1.lc_bas) for i = 1:BS1.nshells]
    lvals2 = [Libcint.CINTcgtos_spheric(i-1, BS2.lc_bas) for i = 1:BS2.nshells]
    Lmax1 = maximum(lvals1)
    Lmax2 = maximum(lvals2)

    # Offset list for each shell, used to map shell index to AO index
    ao_offset1 = [sum(lvals1[1:(i-1)]) for i = 1:BS1.nshells]
    ao_offset2 = [sum(lvals2[1:(i-1)]) for i = 1:BS2.nshells]

    buf = zeros(Cdouble, Lmax1*Lmax2)
    for i in 1:BS1.nshells
        Li = lvals1[i]
        ioff = ao_offset1[i]
        for j in 1:BS2.nshells
            Lj = lvals2[j]
            joff = ao_offset2[j]
            #
            # Call libcint
            libcint_1e!(buf, Cint.([i-1, j+BS1.nshells-1]), lc_atm, natm, lc_bas, nbas, env)
            #
            # Loop through shell block and save unique elements
            for js in 1:Lj
                J = joff + js
                for is in 1:Li
                    I = ioff + is
                    out[I,J] = buf[is + Li*(js-1)]
                end
            end
        end
    end
    return out
end