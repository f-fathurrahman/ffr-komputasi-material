from my_mace import data, modules, tools
import numpy as np
import torch
import torch.nn.functional
from my_e3nn import o3
from matplotlib import pyplot as plt
import ase.io

from ase.visualize import view
from scipy.spatial.transform import Rotation

from my_mace.tools import torch_geometric
torch.set_default_dtype(torch.float64)
import warnings
warnings.filterwarnings("ignore")

# setup some default prameters
z_table = tools.AtomicNumberTable([1, 6, 8])
atomic_energies = np.array([-1.0, -3.0, -5.0], dtype=float)
cutoff = 3

default_model_config = dict(
    num_elements = 3,  # number of chemical elements
    atomic_energies = atomic_energies,  # atomic energies used for normalisation
    avg_num_neighbors = 8,  # avg number of neighbours of the atoms, used for internal normalisation of messages
    atomic_numbers = z_table.zs,  # atomic numbers, used to specify chemical element embeddings of the model
    r_max = cutoff,  # cutoff
    num_bessel = 8,  # number of radial features
    num_polynomial_cutoff = 6,  # smoothness of the radial cutoff
    max_ell = 2,  # expansion order of spherical harmonic adge attributes
    num_interactions = 2,  # number of layers, typically 2
    interaction_cls_first = modules.interaction_classes[
        "RealAgnosticResidualInteractionBlock"
    ],  # interation block of first layer
    interaction_cls = modules.interaction_classes[
        "RealAgnosticResidualInteractionBlock"
    ],  # interaction block of subsequent layers
    hidden_irreps = o3.Irreps("8x0e + 8x1o"),  # 8: number of embedding channels, 0e, 1o is specifying which equivariant messages to use. Here up to L_max=1
    correlation = 3,  # correlation order of the messages (body order - 1)
    MLP_irreps = o3.Irreps("16x0e"),  # number of hidden dimensions of last layer readout MLP
    gate = torch.nn.functional.silu,  # nonlinearity used in last layer readout MLP
)
default_model = modules.MACE(**default_model_config)

# Representing spherical tensors
# a function for Ylms where we evaluate for l=0,1,2.
spherical_harmonics = o3.SphericalHarmonics([0,1,2], True)

# evaulate spherical harmonics on a vector
vector = torch.tensor([1.0, 0.2, 0.75])
print(spherical_harmonics(vector))

# make a list of rotated versions of the vector
rotated_vectors = []
vector = np.array([1.0, 0.2, 0.75])
N = 360
for i in range(N):
    # rotate around the vector [0, 0.7071, 0.7071]
    rotation_matrix = Rotation.from_rotvec(i * 2*np.pi * np.array([0, 0.7071, 0.7071])/360).as_matrix()
    rotated_vectors.append(rotation_matrix @ vector)

# convert to torch tensor
rotated_vectors = torch.tensor(rotated_vectors, dtype=torch.float32)

# compute the spherical harmonics for each vector
spherical_harmonic_values = spherical_harmonics(rotated_vectors)
print('shape of Y_lms array is ', spherical_harmonic_values.shape)

# plot
labels = [[f'l={l}, m={m}' for m in range(-l,l+1)] for l in range(3)]
labels = [x for xs in labels for x in xs] # flatten
plt.plot(spherical_harmonic_values.numpy(), label=labels)
plt.legend()
plt.xlabel('rotation angle')
plt.ylabel('spherical harmonic value')
plt.show()


