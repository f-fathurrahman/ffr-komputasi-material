from pyscf import gto, dft
from pyscf.pbc.tools.pyscf_ase import atoms_from_ase

basis = "STO-3G"
spin_mult = 0
charge = 0

mol = gto.M(
    atom = "N2H4.xyz",
    basis = basis,
    spin = spin_mult,
    charge = charge
)

mf = dft.UKS(mol)
mf.xc = "lda,pw"
mf.max_cycle = 300
mf.conv_tol = 1e-6
mf.verbose = 4
mf.kernel()


