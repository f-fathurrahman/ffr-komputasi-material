#from my_mace import data, modules, tools
import torch
import numpy as np

from ase.build import molecule
from ase.neighborlist import neighbor_list

from my_e3nn import o3
from my_e3nn.o3 import Irreps
#from my_e3nn.nn.models.v2106.gate_points_message_passing import smooth_cutoff


# -----------------------------
# Build example atomic system
# -----------------------------
atoms = molecule("H2O")

positions = torch.tensor(atoms.get_positions(), dtype=torch.float32)
atomic_numbers = atoms.get_atomic_numbers()

cutoff = 3.0


# -----------------------------
# Neighbor list
# -----------------------------
i, j, S = neighbor_list("ijS", atoms, cutoff)

i = torch.tensor(i)
j = torch.tensor(j)

rij = positions[j] - positions[i]
r = rij.norm(dim=1)

rhat = rij / r[:, None]


# -----------------------------
# Radial basis
# -----------------------------
n_radial = 4

def radial_basis(r):
    """Simple Gaussian radial basis."""
    centers = torch.linspace(0, cutoff, n_radial)
    widths = 0.5
    return torch.exp(-(r[:, None] - centers)**2 / widths**2)

radial = radial_basis(r)


# -----------------------------
# Spherical harmonics
# -----------------------------
lmax = 2

Y = o3.spherical_harmonics(
    list(range(lmax + 1)),
    rhat,
    normalize=True,
    normalization="component"
)

# irreps of spherical harmonics
irreps_sh = o3.Irreps.spherical_harmonics(lmax)

print("Spherical harmonic irreps:", irreps_sh)


# -----------------------------
# Initial atom embedding
# -----------------------------
n_atom_features = 16

irreps_atom = Irreps(f"{n_atom_features}x0e")

# simple embedding table
embedding = torch.nn.Embedding(100, n_atom_features)

atom_features = embedding(torch.tensor(atomic_numbers))


# -----------------------------
# Edge features
# -----------------------------
# radial features are scalars (0e)
irreps_radial = Irreps(f"{n_radial}x0e")

# combine radial and angular
tp_radial_sh = o3.FullyConnectedTensorProduct(
    irreps_radial,
    irreps_sh,
    irreps_sh * n_radial
)

edge_features = tp_radial_sh(radial, Y)

print("Edge irreps:", tp_radial_sh.irreps_out)


# -----------------------------
# Message tensor product
# -----------------------------
tp_message = o3.FullyConnectedTensorProduct(
    irreps_atom,
    tp_radial_sh.irreps_out,
    Irreps("16x0e + 16x1o + 8x2e")
)

messages = tp_message(atom_features[j], edge_features)

print("Message irreps:", tp_message.irreps_out)


# -----------------------------
# Aggregate messages
# -----------------------------
n_atoms = len(atoms)

atom_out = torch.zeros(
    n_atoms,
    tp_message.irreps_out.dim
)

for idx in range(len(i)):
    atom_out[i[idx]] += messages[idx]

print("Final atom feature shape:", atom_out.shape)
print("Final atom irreps:", tp_message.irreps_out)