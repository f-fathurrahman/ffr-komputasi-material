# Load Fermi
using Fermi
using Fermi.Integrals
using LinearAlgebra

# ffr: These macro calls will modify some "global" dictionary in Fermi.Options
@molecule {
    O    1.209153654800    1.766411818900   -0.017161397200
    H    2.198480007500    1.797710062700    0.012116171900
    H    0.919788188200    2.458018557000    0.629793883200
}

@set {
    basis sto-3g
    charge 0
    multiplicity 1 # Note that multiplicity must be one for RHF
}

#ffr: The global dict is MyFermi.Options.Current

# It is important to specify `eri_type`. The default will return a Sparse Array.
# The RHF algorithm is faster using this sparse array, but it also gets more complicated
# Here, we are looking for the simplest implementation
aoints = IntegralHelper(eri_type=Chonky())
# ffr: Why not using Dense() ?


# ffr: We only need to pass aoints here
function MyRHF(aoints)

    # Get integrals
    println("Collecting Integrals")
    S = aoints["S"]
    T = aoints["T"]
    V = aoints["V"]
    H = T + V # single particle
    G = 2*aoints["ERI"] - permutedims(aoints["ERI"], (1,3,2,4)) # Coulomb and exchange?
    X = S^(-1/2)

    # Get nuclear repulsion
    Vnuc = aoints.molecule.Vnuc
    
    # Get the number of doubly occupied orbitals
    ndocc = aoints.molecule.Nα
    
    # Get the number of basis functions
    nbf = size(S, 1)
    
    # Create an array for C and set it to zero
    C = zeros(nbf, nbf)
    
    # Get density matrix
    D = C[:,1:ndocc] * (C[:,1:ndocc])'
    
    # Starts iterations!
    ΔE = 1.0 # arbitrary, just to start the loop
    Eold = 0.0
    Enew = 0.0
    ϵ = zeros(nbf)
    
    println("Starting Iterations!")
    
    while ΔE > 1e-8 

        Eold = Enew
        
        # Build Fock matrix
        F = similar(H)
        F .= H # Don't do F = H !!
        
        for μ in 1:nbf, ν in 1:nbf, ρ in 1:nbf, σ in 1:nbf
            F[μ,ν] += G[μ,ν,ρ,σ]*D[σ,ρ]
        end
        
        # Tarsnform F
        tF = X'*F*X
        
        # Diagonalize F
        ϵ, tC = LinearAlgebra.eigen(Symmetric(tF), sortby=x->x)
        
        # Backtransform C
        C = X*tC
        
        # Update density matrix
        D = C[:,1:ndocc] * (C[:,1:ndocc])'
        
        # Compute energy
        Enew = Vnuc
        for μ in 1:nbf, ν in 1:nbf
            # Watch out! This portion cannot be all the way inside the loop
            Enew += 2*H[μ,ν]*D[μ,ν]
            for ρ in 1:nbf, σ in 1:nbf
                Enew += G[μ,ν,ρ,σ]*D[σ,ρ]*D[μ,ν]
            end
        end
              
        # Compute ΔE
        ΔE = abs(Enew - Eold)
        
        # Print some msg!
        println("New energy: $Enew  - ΔE = $ΔE")
    end
    
    # Return energy and orbitals
    return Enew, ϵ, C
end