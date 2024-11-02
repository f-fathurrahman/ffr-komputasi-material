import numpy as np
import jax
from jax_dft import utils, scf
import matplotlib.pyplot as plt

jax.config.update("jax_enable_x64", True)
print(f"JAX devices: {jax.devices()}")

num_electrons = 2
grids = np.arange(-256, 257)*0.08
external_potential = utils.get_atomic_chain_potential(
    grids=grids,
    locations=np.array([-0.8, 0.8]), # needs to be within the grids
    nuclear_charges=np.array([1.0, 1.0]),
    interaction_fn=utils.exponential_coulomb
)

density, total_eigen_energies, _ = scf.solve_noninteracting_system(
    external_potential,
    num_electrons=num_electrons, grids=grids
)

print(f"total energy: {total_eigen_energies}")

plt.plot(grids, density, label="density")
plt.plot(grids, external_potential, label="potential")
plt.legend(loc=0)
plt.grid(True)
plt.savefig("IMG_non_interacting.png", dpi=150)
plt.show()
