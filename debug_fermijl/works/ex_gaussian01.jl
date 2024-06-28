using MyGaussianBasis

# %% Construct a basis set with a given molecule
bs_HCl = BasisSet("sto-3g", """
H         0.00      0.00     0.00
Cl        0.76      0.00     0.00""")

# The actual basis functions for each atom
bs_HCl.basis[1] # for H of each of type Vector{BasisFunction}
bs_HCl.basis[2] # for Cl

Nbasis_H = length(bs_HCl.basis[1])
Nbasis_Cl = length(bs_HCl.basis[2])

for (i,b) in enumerate(bs_HCl.basis[1])
    println("\nBasis number: $i")
    display(b)
end

for (i,b) in enumerate(bs_HCl.basis[2])
    println("\nBasis number: $i")
    display(b)
end