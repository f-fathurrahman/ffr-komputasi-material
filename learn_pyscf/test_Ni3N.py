from pyscf import gto, dft, scf

spin_mult = 1
charge = 0
mol = gto.M(
    atom = "Ni3N.xyz",
    spin = spin_mult,
    charge = charge
)
mol.basis = "gth-dzvp-molopt-sr"
mol.pseudo = "gth-pbe"
mol.build()

mf = dft.UKS(mol)
mf = scf.addons.smearing_(mf, sigma=0.1, method='fermi')
mf.xc = "pbe"
mf.max_cycle = 300
mf.conv_tol = 1e-8
mf.verbose = 4 # optimal info
mf.kernel()

from pyscf.geomopt.geometric_solver import optimize as geometric_opt
optim_geo1 = geometric_opt(mf)



