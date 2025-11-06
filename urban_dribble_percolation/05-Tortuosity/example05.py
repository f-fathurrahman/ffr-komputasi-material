# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3.8.8 ('base')
#     language: python
#     name: python3
# ---

# %% [markdown]
# #### Install `Dribble` and its requirements

# %%
try:
    import dribble
except ImportError:
    # %pip install git+https://github.com/atomisticnet/dribble.git

# %% [markdown]
# # Tortuosity
#
# The *tortuosity* is a measure of the *detour* that a percolating diffusion pathway takes.  If the length of the diffusion pathway between two periodic images of the same site is $L$ and the total distance between the periodic images is $D$, the tortuosity is
#
# \begin{equation}
#   \tau = \frac{L}{D}
#   \quad .
# \end{equation}
#
# This means, a tortuosity of $1$ is ideal (no detour), and the larger the tortuosity becomes, the greater the detour that has to be taken.
#
# ## Input File
#
# Here, we will determine the tortuosity for the nearest-neighbor site percolation problem of [Example 0](../00-Command-Line-Usage/example00.ipynb). Once again, the same general input file can be used:

# %%
# %%writefile input-bond-rule.json
{
    "structure": "Cu-fcc.vasp",
    "formula_units": 1.0,
    "sublattices": {
        "A": {
            "description": "Copper sites",
            "sites": {"species": ["Cu"]},
            "initial_occupancy": {"Vac": 1.0}
        }
    },
    "bonds": [
        {
            "sublattices": ["A", "A"],
            "bond_rules": [["NearestNeighborBR"]]
        }
    ],
    "percolating_species": ["Cu"],
    "flip_sequence": [["Vac", "Cu"]]
}

# %%
# %%writefile Cu-fcc.vasp
FCC Structure
3.6
     0.0      0.5      0.5
     0.5      0.0      0.5
     0.5      0.5      0.0
Cu
1
direct
     0.0      0.0      0.0 Cu

# %% [markdown]
# See [Example 0](../00-Command-Line-Usage/example00.ipynb) for a detailed discussion of the input file.
#
# ## Calculation using the Command Line Tool
#
# The tortuosity is computed by `Dribble`'s command line tool if the `-t` (or `--tortuosity`) flag is present.  Once again, care has to be taken to select a sufficiently large supercell.  Note, however, that computing the tortuosity can be extremely time consuming, as it involves determining the shortest pathway between two endpoints on a graph (the network of percolating pathways).  Therefore, we will use a much smaller supercell and fewer MC samples for this example:

# %%
# ! dribble input-bond-rule.json --supercell 4 4 4 -t --samples 10

# %% [markdown]
# The results of the tortuosity simulation are written to a file named `percol.tortuosity`.  Let's take a look at the generated output:

# %%
# !head -n20 percol.tortuosity

# %% [markdown]
# The first column in the output file is the number of occupied sites, the second column is the fraction of occupied sites relative to the total number of sites, and the third column is the corresponding tortuosity.  **Below the percolation threshold ($x_c\approx{}0.199$) the tortuosity is not defined, so the values are not meaningful.**  For concentrations that are never found to be percolating, the tortuosity is simply set to infinity ("inf").  However, in the implemented Monte Carlo method structures may sometimes "by coincidence" become percolating even below the percolation threshold, and these values should be ignored.
#
# Let's plot the tortuosity for $x>x_c$:

# %%
# %matplotlib inline
from matplotlib import pyplot as plt
import numpy as np
data = np.loadtxt("percol.tortuosity")
fig, ax = plt.subplots(figsize=(8,5))

x_c = 0.199

ax.set_xlim(0.199,1)
ax.set_xlabel("Occupancy ($p$)", fontsize=16)
ax.set_ylabel(r"Tortuosity (\tau)", fontsize=16, color="tab:blue")
ax.tick_params(labelsize=16)
ax.tick_params(axis='y', labelcolor="tab:blue")

ax.plot(data[:,1], data[:,2], linewidth=3, color="tab:blue")
plt.show()

# %% [markdown]
# As seen in the plot, tortuosity is not yet very well converged as we only used 10 Monte Carlo samples.  Additionally, the values near the percolation threshold are likely not yet converged with the cell size.  However, the overall trend is clear: For concentrations just above the percolation threshold the tortuosity is $\tau>1.6$, i.e., diffusion pathways lead to detours of more than 60% of the diffusion distance.  Only above site concentrations of $\sim{}0.7$ the tortuosity approaches $1.0$ meaning that the diffusion pathways span the structure without detour.
