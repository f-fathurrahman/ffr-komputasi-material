import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import ase.io
atoms_list = ase.io.read("DATASET_OTHERS/TiAl_gabung.xyz@:")
atoms = atoms_list[0]

positions = atoms.get_positions()
# Create data for the plot
x = positions[:,0]
y = positions[:,1]
z = positions[:,2]
Natoms = len(atoms)
sizes = np.random.rand(Natoms) * 300  # Generate random sizes for each point

# Create the figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create the scatter plot, adjusting point size with 's'
scatter = ax.scatter(x, y, z, s=sizes)

# Set labels for the axes
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')

# Show the plot
plt.show()