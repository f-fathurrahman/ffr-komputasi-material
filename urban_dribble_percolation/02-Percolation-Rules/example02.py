import sys
sys.path.append("../")

# # Bond and Site Rules
#
# ``Dribble`` implements a number of different **bond rules** and **site rules** that determine the
# criteria under which sites bond with each other.  In general, these rules differ from material to
# material and from application to application.  
#
# For example, Li conduction in cation-disordered Li transition-metal (TM) oxides with rocksalt
# structure is known to proceed via so-called **0-TM** diffusion channels [1].  The rocksalt structure
# consist of two FCC sublattices, one for the cations (Li and TM in this case) and one for the oxygen atoms.
# Since the oxygen sites are all equivalent, we are again interested in percolation on the FCC lattice as in [example 00](../00-Command-Line-Usage/example00.ipynb).  However, instead of just a nearest-neighbor bond rule, we now have to encode an appropriate rule for 0-TM channels.
#
# [1] [Lee et al., Science 31, 2014, 519-522](http://doi.org/10.1126/science.1246432)


# ## Bond Rules
#
# One way to think of 0-TM channels is in terms of a *bond* criterion between two Li sites:
#
# A 0-TM channel between two neighboring Li sites $i$ and $j$ exists if 
#
# 1. At least 2 sites $k$ and $l$ that are nearest neighbors of both $i$ and $j$ are also Li sites; and
# 2. $k$ and $l$ are themselves nearest neighbors.
#
# We can express this condition using the `MinCommonNNNeighborsBR` bond rule provided by `dribble` (assuming that the cation sublattice is named "cations"):
#
# ```
#     "bonds": [
#         {
#             "sublattices": ["cations", "cations"],
#             "bond_rules": [["MinCommonNNNeighborsBR", {"num_neighbors": 2}]]
#         }
#     ],
# ```
#
# The entire `dribble` input file with this bond rule is as follows:

# Note how the oxygen sites are set to be ignored in the above input file.
#
# ## Site Rules
#
# Alternatively to the above *bond rule* definition of 0-TM channels, we can also think of 0-TM channels as tetrahedral sites that are coordinated at all four faces by Li sites.  The *normal* cation sites in the rocksalt structure are *octahedral* sites, but the Li diffusion takes place via a *tetrahedral* intermediate.  To use such a **site rule** criterion, we require a structure file that also includes the tetrahedral sites (corresponding to the calcium fluorite structure).  Then we can define tetrahedral sites to be only accessible when coordinated by 4 Li sites.  The *sublattice* block for the tetrahedral sublattice looks as follows:
#
# ```
#         "tet": {
#             "description": "Tetrahedral site",
#             "sites": [3, 4, 5, 6],
#             "initial_occupancy": {"Vac": 1.0},
#             "site_rules": [
#                 ["NeighborShellSR",
#                  {"stable_nb_shells": [[
#                      {"oct": [{"min": 4, "species": ["Li"]}]}
#                    ]]
#                  }
#                 ]
#             ]
#         },
#
# ```
#
# There are a few differences to the previous input file:
#
# - The tetrahedral sites are defined explicitly by their occurance in the structure file (sites 3-6) instead of by species;
# - The **site rule** is defined in a new sub-block using neighbor shell occupations (implemented as `NeighborShellSR` site rule in `dribble`).  A neighbor shell site rule defines under which neighbor shell conditions the sites of a given sublattice become accessible.  Here, we require the first neighbor shell of each tetrahedral site to contain at least (`min`) 4 octahedral (`oct` sublattice) Li sites.
#
# The complete `dribble` input file is as follows:


writefile input-site-rule.json
{
    "structure": "LiMO2+tet-rotated.vasp",
    "formula_units": 1,
    "cutoff": 2.0,
    "sublattices": {
        "tet": {
            "description": "tetrahedral site",
            "sites": [3, 4, 5, 6],
            "initial_occupancy": {"Vac": 1.0},
            "site_rules": [
                ["NeighborShellSR",
                 {"stable_nb_shells": [[
                     {"oct": [{"min": 4, "species": ["Li"]}]}
                   ]]
                 }
                ]
            ]
        },
        "oct": {
            "description": "octahedral site",
            "sites": [1, 2],
            "initial_occupancy": {"TM": 1.0}
        },
        "oxygen": {
            "description": "oxygen sites",
            "sites": {"species": ["O"]},
            "ignore": true
        }
    },
    "bonds": [{"sublattices": ["oct", "tet"]}],
    "percolating_species": ["Li", "Vac"],
    "static_species": ["Vac"],
    "flip_sequence": [["TM", "Li"]]
}


# Note three more differences compared to the *bond rule* input file:
#
# - An explicit **cutoff** of 2.0 Å is specified for the neighbor shell detection;
# - **Bonds** are now between octahedral (`oct`) and tetrahedral (`tet`) sites; and
# - Tetrahedral sites are occupied by species Vac which is defined to be a **static species**, i.e., the Vac sites do not change during the percolation simulation.
#
# ## Percolation Simulations
#
# Now let's compare the results obtained using the two different percolation rules.
#
# 1. Using the *bond rules*:


# ! dribble input-bond-rule.json --supercell 6 6 6 -p


# 2. Using the *site rules*:


# ! dribble input-site-rule.json --supercell 6 6 6 -p


# ## Discussion
#
# - The literature value for the percolation threshold is $\approx 1.09$ Li per formula unit, and both approaches get this about right.
# - The computational cost of the *site rule* implementation is slightly higher as the simulation involves a greater total number of sites.  However, some diffusion channels might be more naturally defined in terms of intermediate sites than in terms of bonds.


