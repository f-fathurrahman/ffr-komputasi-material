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
# # Wrapping Probability
#
# One way to estimate the percolation threshold is by evaluating the probability of *periodically wrapping* diffusion pathways existing as a function of the concentration of *occupied sites* (see also [Example 0](../00-Command-Line-Usage/example00.ipynb)).  A diffusion pathway is periodically wrapping if it connects a site with one of its periodic images.
#
# ## Input File
#
# Here, we will determine the periodic wrapping probability for the nearest-neighbor site percolation problem of [Example 0](../00-Command-Line-Usage/example00.ipynb). The same general input file can be used:

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
# The periodic wrapping probability is computed by `Dribble`'s command line tool if the `-w` (or `--pwrap`) flag is present.  Once again, care has to be taken to select a sufficiently large supercell:

# %%
# ! dribble input-bond-rule.json --supercell 10 10 10 -w

# %% [markdown]
# Per default, `Dribble` returns a binomial convolution of the Monte Carlo results which helps to improve convergence with the number of Monte Carlo samples (500 per default).  These results are written to a file named `percol.wrap`.  The raw data can be requested with the flag `--save-raw`.
#
# Let's take a look at the generated output file:

# %%
# !head percol.wrap

# %% [markdown]
# The first column in the output file contains the relative site occupancy, i.e., the fraction of sites that is occupied.  The second column gives the probability that periodic wrapping is first detected at the given site occupancy, and the third column is the cumulative, i.e., the probability that periodic wrapping occurred at any concentration up to a given one (the integral of the second column).
#
# Maybe a plot makes this easier to grasp:

# %%
# %matplotlib inline
from matplotlib import pyplot as plt
import numpy as np
data = np.loadtxt("percol.wrap")
fig, ax = plt.subplots(figsize=(8,5))
ax.set_xlim(0,1)
ax.set_xlabel("Occupancy ($p$)", fontsize=16)
ax.set_ylabel(r"Wrapping Probability ($P_{\rm{wrap}}$)", fontsize=16, color="tab:blue")
ax.tick_params(labelsize=16)
ax.tick_params(axis='y', labelcolor="tab:blue")

ax2 = ax.twinx()
ax2.set_ylabel(r"Cumulative Wrapping Probability", fontsize=16, color="tab:red")
ax2.tick_params(labelsize=16, labelcolor="tab:red")

ax.plot(data[:,0], data[:,1], linewidth=3, color="tab:blue")
ax2.plot(data[:,0], data[:,2], linewidth=3, color="tab:red")
plt.show()

# %% [markdown]
# As seen in the plot, the periodic wrapping probability has a peak at the percolation threshold $x_c\approx{}0.199$.  With increasing cell size, the cumulative approaches a step function that is $0$ below the percolation threshold and $1$ above.
