# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# ## Setup

# %%
include("setup_works.jl")

# %% [markdown]
# ## Load library

# %%
# Load Fermi
using MyFermi
using MyFermi.Integrals

# %%
using LinearAlgebra

# %% [markdown]
# These macro calls will modify some "global" dictionary in `MyFermi.Options`:

# %%
MyFermi.Options.Current

# %%
@molecule {
    O    1.209153654800    1.766411818900   -0.017161397200
    H    2.198480007500    1.797710062700    0.012116171900
    H    0.919788188200    2.458018557000    0.629793883200
}

# %% [markdown]
# After calling `@molecule` macro, the dict `MyFermi.Options.Current` will be modified:

# %%
MyFermi.Options.Current

# %% [markdown]
# Setting basis set:

# %%
@set {
    basis sto-3g
    charge 0
    multiplicity 1 # Note that multiplicity must be one for RHF
}

# %%
MyFermi.Options.Current

# %% [markdown]
# ## IntegralHelper

# %% [markdown]
# IntegralHelper seems to hold all important information to carry out RHF calculations.

# %% [markdown]
# It is important to specify `eri_type`. The default will return a `SparseArray`.
# The RHF algorithm is faster using this sparse array, but it also gets more complicated
# Here, we are looking for the simplest implementation

# %%
aoints = IntegralHelper(eri_type=Chonky())

# %% [markdown]
# The integrals are store in dict `cache`. Initially, cache is empty:

# %% editable=true slideshow={"slide_type": ""}
keys(aoints.cache)

# %% [markdown]
# Accessing `aoints` like a dict will call various wrappers to compute molecular integrals. As an example, we will access overlap integrals:

# %%
S = aoints["S"]

# %% [markdown]
# When getindex of IntegralHelper is called for the **first time**, it will call compute! which wraps calls to compute_S!, compute_T!, etc and save the result in `aoints.cache`:

# %%
keys(aoints.cache)

# %% [markdown]
# Functions to calculate molecular integrals:
# - compute_S!
# - compute_T!
# - compute_V!
# - compute_ERI!

# %% [markdown]
# Now, we will call other molecular integrals:

# %%
aoints["T"] # kinetic

# %%
aoints["V"] # nuclear

# %%
keys(aoints.cache)

# %%
aoints["S"]

# %% [markdown]
# If we call various `Integrals.compute`s functions directly, it will force recalculation and the results are also cached.

# %%
MyFermi.Integrals.compute_ERI!(aoints);

# %%
keys(aoints.cache)

# %%
aoints["ERI"] |> size

# %% [markdown]
# Molecule information is also stored in `IntegralHelper`:

# %%
aoints.molecule

# %%
aoints |> typeof |> fieldnames

# %%
aoints.molecule |> typeof

# %%
using MyMolecules

# %%
MyMolecules.nuclear_repulsion(aoints.molecule.atoms)

# %% [markdown]
# ## MyRHF Implementation

# %% editable=true slideshow={"slide_type": ""}
# ffr: We only need to pass aoints here
function MyRHF_v1(aoints)

    # Get integrals
    println("Collecting Integrals")
    S = aoints["S"]
    T = aoints["T"]
    V = aoints["V"]
    H = T + V # single particle
    G = 2*aoints["ERI"] - permutedims(aoints["ERI"], (1,3,2,4)) # Coulomb and exchange?
    X = S^(-1/2)

    # Get nuclear repulsion
    #Vnuc = aoints.molecule.Vnuc
    Vnuc = MyMolecules.nuclear_repulsion(aoints.molecule.atoms)
    
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

    F = similar(H) # ffr: pull this out of the loop?
    
    println("Starting Iterations!")
    
    while ΔE > 1e-8 

        Eold = Enew
        
        # Build Fock matrix
        #F = similar(H) # ffr: pull this out of the loop?
        F .= H # Don't do F = H !!
        
        for μ in 1:nbf, ν in 1:nbf, ρ in 1:nbf, σ in 1:nbf
            F[μ,ν] += G[μ,ν,ρ,σ]*D[σ,ρ]
        end
        
        # Transform F
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

# %%
@time res1 = MyRHF_v1(aoints)

# %%
@time res1 = MyRHF_v1(aoints)

# %%
@time res1 = MyRHF_v1(aoints)

# %%

# %%
