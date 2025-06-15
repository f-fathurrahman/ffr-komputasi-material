# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Preparation

# %%
from ase.build import bulk
from my_emt import MyEMT
from my_phonons import MyPhonons
# Setup crystal and EMT calculator
atoms = bulk("Al", "fcc", a=4.05)
# Phonon calculator
N = 7
phonon = MyPhonons(atoms, MyEMT(), supercell=(N, N, N), delta=0.05)
phonon.run()

# %% [markdown]
# ## Read forces, calculate phonons

# %%
import numpy as np

# %% [markdown]
# Some hardcoded parameters:

# %%
method = "Frederiksen"
symmetrize = 3
acoustic = True
cutoff = None
born = False

# %%
# Number of atoms
natoms = len(phonon.indices)
print("natoms = ", natoms)

# Number of unit cells
N = np.prod(phonon.supercell)
print("Number of unit cells = ", N)

# %% [markdown]
# Number of unit cells comes from here:

# %%
7**3

# %%
# Matrix of force constants as a function of unit cell index in units
# of eV / Ang**2
C_xNav = np.empty((natoms*3, N, natoms, 3), dtype=float)

# Loop over all atomic displacements and calculate force constants
for i, a in enumerate(phonon.indices): # loop for all atoms
    for j, v in enumerate("xyz"): # displacements
        #
        # Atomic forces for a displacement of atom a in direction v
        basename = "%d%s" % (a, v)
        fminus_av = phonon.cache[basename + "-"]["forces"]
        fplus_av = phonon.cache[basename + "+"]["forces"]

        if method == "frederiksen":
            fminus_av[a] -= fminus_av.sum(0)
            fplus_av[a] -= fplus_av.sum(0)

        # Finite difference derivative (central difference for 1st derivative)
        C_av = fminus_av - fplus_av
        C_av /= 2 * phonon.delta
        print("shape of C_av = ", C_av.shape)
        #
        # Slice out included atoms
        C_Nav = C_av.reshape((N, len(phonon.atoms), 3))[:, phonon.indices]
        index = 3*i + j
        C_xNav[index] = C_Nav

# %% [markdown]
# The results are in `C_xNav`.

# %%
C_xNav[0,3,0,:] # (3, 7*7*7, Natoms, 3) ?

# %%
# Make unitcell index the first and reshape
C_N = C_xNav.swapaxes(0, 1).reshape((N,) + (3 * natoms, 3 * natoms))

# %%
C_N.shape

# %%
C_N[0]

# %%
from m_debug_phonon import ph_symmetrize, ph_acoustic
# Symmetrize force constants
print("symmetrize = ", symmetrize)
if symmetrize:
    for i in range(symmetrize):
        # Symmetrize
        C_N = ph_symmetrize(phonon, C_N)
        # Restore acoustic sum-rule
        if acoustic:
            ph_acoustic(phonon, C_N)
            #phonon.acoustic(C_N)
        else:
            break

# %%
C_N[0]

# %%
# Store force constants and dynamical matrix
phonon.C_N = C_N
phonon.D_N = C_N.copy()

# Add mass prefactor
m_a = phonon.atoms.get_masses()
phonon.m_inv_x = np.repeat(m_a[phonon.indices]**-0.5, 3)
M_inv = np.outer(phonon.m_inv_x, phonon.m_inv_x)
for D in phonon.D_N:
    D *= M_inv
