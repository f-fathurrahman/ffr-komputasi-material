from pyscf import gto, dft, scf, df

mol = gto.M(
    atom = "N2H4.xyz",
    spin = 0,
    charge = 0
)
mol.basis = "sto-3g"
mol.build()

mf = dft.UKS(mol)
mf.xc = "pbe"
mf.max_cycle = 300
mf.verbose = 4 # optimal info

mf = mf.density_fit()
mf.with_df = df.DF(mol).build()
mf.kernel()

from pyscf.geomopt.geometric_solver import optimize as geometric_opt
optim_geo1 = geometric_opt(mf)



