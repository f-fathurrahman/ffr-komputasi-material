"""
Test case for Z92 Schrodinger mode against analytic hydrogenic levels.

This test validates the `xc_functional='Schrodinger'` mode where both XC and
Hartree contributions are disabled in the Hamiltonian build path.
"""

import time
import numpy as np

from atom.solver import AtomicDFTSolver

def hydrogenic_analytic_energy(n: int, z: int) -> float:
    return -(z**2) / (2.0 * (n**2))

def z92_schrodinger_occupied_abs_errors(results, z: int = 92) -> np.ndarray:
    """Per-occupied-state |E_numeric - E_analytic| (Ha) for hydrogenic comparison."""
    occ_info = results["occupation_info"]
    n_occ = occ_info.n_states
    full_eigs = np.asarray(results["full_eigen_energies"], dtype=float)
    if full_eigs.shape[0] < n_occ:
        raise ValueError(
            f"full_eigen_energies length {full_eigs.shape[0]} < n occupied states {n_occ}"
        )

    occ_n = np.asarray(occ_info.occ_n, dtype=int)
    numeric = full_eigs[:n_occ]
    analytic = np.array(
        [hydrogenic_analytic_energy(n=int(occ_n[i]), z=z) for i in range(n_occ)],
        dtype=float,
    )
    return np.abs(numeric - analytic)


print("Running Z92 Schrodinger calculation...")
start_time = time.time()
solver = AtomicDFTSolver(
    atomic_number           = 92,
    n_electrons             = 92,
    xc_functional           = "Schrodinger",
    domain_size             = 40.0,
    finite_element_number   = 12,
    polynomial_order        = 20,
    quadrature_point_number = 60,
    mesh_type               = "exponential",
    mesh_concentration      = 101.0,
    scf_tolerance           = 1e-20,
    verbose                 = True,
    all_electron_flag       = True,
    use_oep                 = False,
    use_preconditioner      = True,
)
results = solver.solve(save_full_spectrum=True, use_warm_start=False)
elapsed = time.time() - start_time
print(f"Computation time: {elapsed:.2f} seconds")

z = 92

occ_info = results["occupation_info"]
n_occ    = occ_info.n_states
full_eigs = np.asarray(results["full_eigen_energies"], dtype=float)
if full_eigs.shape[0] < n_occ:
    raise ValueError(
        f"full_eigen_energies length {full_eigs.shape[0]} < n occupied states {n_occ}"
    )

occ_n = np.asarray(occ_info.occ_n, dtype=int)
occ_l = np.asarray(occ_info.occ_l, dtype=int)
f_occ = np.asarray(occ_info.occ_spin_up_plus_spin_down, dtype=float)
degen = 2 * (2 * occ_l + 1)

abs_err = z92_schrodinger_occupied_abs_errors(results, z=z)
numeric = full_eigs[:n_occ]
analytic = np.array(
    [hydrogenic_analytic_energy(n=int(occ_n[i]), z=z) for i in range(n_occ)],
    dtype=float,
)
max_err = float(np.max(abs_err))

print("\nPer-state comparison (all occupied KS states, hydrogenic E_n = -Z^2/(2n^2)):")
print(
    f"{'#':<4} {'n':<4} {'l':<4} {'occ':<6} {'g':<4} "
    f"{'numeric (Ha)':<22} {'analytic (Ha)':<22} {'|err|':<12}"
)
print("-" * 92)
for i in range(n_occ):
    print(
        f"{i:<4} {occ_n[i]:<4} {occ_l[i]:<4} {f_occ[i]:<6.4g} {degen[i]:<4} "
        f"{numeric[i]:<22.13f} {analytic[i]:<22.13f} {abs_err[i]:<12.3e}"
    )
print(f"\nMax abs error (over {n_occ} occupied states): {max_err:.6e} Ha")

tolerance = 1e-5
if max_err < tolerance:
    print(f"Z92 Schrodinger eigenvalues (max err < {tolerance:.1e} Ha)")
