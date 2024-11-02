import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import jax
import jax.numpy as jnp
from jax_dft import scf, utils

# Set the default dtype as float64
jax.config.update("jax_enable_x64", True)

COLORS = [
    "#0072b2",
    "#de8f05",
    "#009e73",
    "#cc79a7",
    "#a24f00",
    "#9467bd",
    "#56b4e9",
    "#bcbd22",
    "#7f7f7f",
]

def set_matplotlib_style():
    """Sets the matplotlib style for the colab notebook."""
    mpl.rcParams["image.cmap"] = "inferno"
    # Set width and size for lines and markers.
    mpl.rcParams["lines.linewidth"] = 2.5
    mpl.rcParams["lines.markersize"] = 9
    mpl.rcParams["lines.markeredgewidth"] = 0
    # Set fontsize.
    mpl.rcParams["font.size"] = 18
    mpl.rcParams["axes.labelsize"] = 20
    mpl.rcParams["axes.titlesize"] = 20
    mpl.rcParams["axes.formatter.useoffset"] = False
    mpl.rcParams["legend.fontsize"] = 14
    mpl.rcParams["xtick.labelsize"] = 14
    mpl.rcParams["ytick.labelsize"] = 14
    # default plot colors recommended by @zan
    # adapted from the seaborn colorblind colors.
    # https://seaborn.pydata.org/tutorial/color_palettes.html
    # Purposely avoids red both for colorblind reasons
    # and because red often reads as more important/bad or more eye-catching than
    # other colors.
    # These colors have been vetted by an individual with red/green color
    # confusion.
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=COLORS)
    mpl.rcParams["savefig.dpi"] = 120
    mpl.rcParams["savefig.bbox"] = "tight"

set_matplotlib_style()


def show_density_potential(
    grids, density, potential, do_show=True, grey=False, axs=None):
    if axs is None:
        _, axs = plt.subplots(nrows=2)
    axs[0].plot(grids, density, c="0.5" if grey else COLORS[0])
    axs[1].plot(grids, potential, c="0.5" if grey else COLORS[1])
    axs[0].set_ylabel(r"$n(x)$")
    axs[1].set_ylabel(r"$v(x)$")
    axs[1].set_xlabel(r"$x$")
    if do_show:
        plt.show()


grids = np.linspace(-5, 5, 201)
dx = utils.get_dx(grids)

# Quantum Harmonic Oscillator: V(x) = 0.5 * k * x**2
# k = 1
# The noninteracting ground state energy is 0.5 Hartree.

qho_potential = 0.5 * grids ** 2
# solve for one electron
qho_density, qho_energy, _ = scf.solve_noninteracting_system(
    qho_potential,
    num_electrons=1,
    grids=grids
)

print(f"total energy: {qho_energy}")
show_density_potential(grids, qho_density, qho_potential, grey=True)


# Perturbed QHO
perturbed_potential = qho_potential + np.exp(-(grids - 0.5) ** 2 / 0.04)
perturbed_density, perturbed_energy, _ = scf.solve_noninteracting_system(
    perturbed_potential,
    num_electrons=1,
    grids=grids
)

print(f"total energy: {perturbed_energy}")
_, axs = plt.subplots(nrows=2)
show_density_potential(
    grids, qho_density, qho_potential, grey=True, do_show=False, axs=axs)
show_density_potential(
    grids, perturbed_density, perturbed_potential, axs=axs)


#
# Adjust potential from loss
#
# Note the use of `jnp` not `np` here.
def density_loss(output, target):
    return jnp.sum((output - target) ** 2) * dx

def energy_loss(output, target):
    return (output - target) ** 2

print(f"Current density loss {density_loss(perturbed_density, qho_density)}")
print(f"Current energy loss {energy_loss(perturbed_energy, qho_energy)}")
print(f"Current total loss {density_loss(perturbed_density, qho_density) + energy_loss(perturbed_energy, qho_energy)}")


# Loss function, given a potential
def loss_fn(potential):
    # get density and energy
    density, energy, _ = scf.solve_noninteracting_system(
        potential,
        num_electrons=1,
        grids=grids
    )
    # compute density loss and energy loss (difference from reference energy and density)
    return density_loss(density, qho_density) + energy_loss(energy, qho_energy)

# Gradient of loss function
grad_fn = jax.jit(jax.grad(loss_fn))  # Compile with jit for fast grad.


plt.plot(grids, grad_fn(perturbed_potential), "--", c=COLORS[2])
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{\partial L_n}{\partial v}$")
plt.show()

# Now we have the gradient. Let"s update the potential from
# the gradient of loss with respect to the potential.
# v   <--   v - epsilon * dL/dv

potential = perturbed_potential
loss_history = []
potential_history = []
record_interval = 1000
for i in range(5001):
    if i % record_interval == 0:
        loss_value = loss_fn(potential)
        print(f"step {i}, loss {loss_value}")
        loss_history.append(loss_value)
        potential_history.append(potential)
    potential -=  30 * grad_fn(potential)

history_size = len(loss_history)

plt.plot(np.arange(history_size) * record_interval, loss_history)
plt.axhline(y=0, color="0.5", ls="--")
plt.xlabel("step")
plt.ylabel(r"$L$")
plt.show()

_, axs = plt.subplots(
    nrows=2, ncols=history_size, figsize=(2.5 * history_size, 4),
    sharex=True, sharey="row"
)

for i, ax in enumerate(axs[0]):
    ax.plot(grids, qho_density, c="0.5")
    density, _, _ = scf.solve_noninteracting_system(
        potential_history[i],
        num_electrons=1,
        grids=grids
    )
    ax.plot(grids, density, "--", c=COLORS[0])
    ax.set_title(rf"$L=${loss_fn(potential_history[i]):1.1e}")

for i, ax in enumerate(axs[1]):
    ax.plot(grids, qho_potential, c="0.5")
    ax.plot(grids, potential_history[i], "--", c=COLORS[1])

# Zoom in the potential.
axs[1][0].set_xlim(-2, 2)
axs[1][0].set_ylim(0.01, 3)
axs[0][0].set_ylabel(r"$n(x)$")
axs[1][0].set_ylabel(r"$v(x)$")
plt.show()

optimized_potential = potential_history[-1]
optimized_density, optimized_total_eigen_energies, _ = (
    scf.solve_noninteracting_system(
        optimized_potential,
        num_electrons=1,
        grids=grids))


print(f"total energy: {optimized_total_eigen_energies}")

_, axs = plt.subplots(nrows=2)
axs[0].plot(grids, optimized_density - qho_density, c=COLORS[0])
axs[0].set_ylabel(r"$\Delta n(x)$")
axs[1].plot(grids, optimized_potential - qho_potential, c=COLORS[1])
axs[1].set_ylabel(r"$\Delta v(x)$")
axs[1].set_xlabel(r"$x$")
plt.show()

_, axs = plt.subplots(nrows=2)
show_density_potential(
    grids, qho_density, qho_potential, grey=True, do_show=False, axs=axs)
show_density_potential(
    grids, optimized_density, optimized_potential, axs=axs)