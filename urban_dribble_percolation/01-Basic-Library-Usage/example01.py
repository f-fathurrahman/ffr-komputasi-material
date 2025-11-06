import sys
sys.path.append("../")

# 
# Basic Usage of the ``Dribble`` Library
#
# For most applications, the three object classes `Input`, `Percolator`, and `Lattice` are needed.
from dribble import Input, Lattice, Percolator

# 
# The `Input` class is used to load input fles in the ``JSON`` format and takes care of setting up sublattices and percolation rules.  An example ``Dribble`` input file named ``input-bond-rule.json`` is in the present directory. 
#
# Load the input file with: 
inp = Input.from_file("input-bond-rule.json")

# 
# The input file also contains a path to a structure file, so that the lattice for the
# simulation can be constructed just based on the Input object.  The results of
# the Monte-Carlo simulation will depend on the size of the simulation cell, so that
# we here generate a 4x4x4 supercell:
lat = Lattice.from_input_object(inp, supercell=[4, 4, 4])

# 
# Next, we set up an instance of the `Percolator` class that takes care of actual
# Monte Carlo simulations.  It only requires an `Input` object and a `Lattice` object as input.
percolator = Percolator.from_input_object(inp, lat, verbose=True)

# 
# Now, we can use the `percolator` to compute various quantities.
# The most basic example is the percolation threshold, which can be obtained with
# the following command:

# 
(pc_site_any, pc_site_two, pc_site_all, pc_bond_any, pc_bond_two, pc_bond_all
) = percolator.percolation_point(inp.flip_sequence, samples=100)

# 
# If everything went well, the Monte-Carlo simulation should give as a result an
# average percolating composition of approximately Li<sub>1.1</sub>TM<sub>0.9</sub>.
# Note that the exact result may vary a little due to finite size effects and the limited
# number of samples (here only 100) in the simulation.  Also note that the composition does
# not contain any oxygen, as the oxygen sublattice was defined to be static in the input file.
