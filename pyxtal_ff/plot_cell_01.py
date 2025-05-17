import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Define lattice vectors for a simple cubic unit cell
a = np.array([1, 0, 0])
b = np.array([0, 1, 0])
c = np.array([0, 0, 1])

# Calculate unit cell vertices
vertices = [
    np.array([0, 0, 0]), a, b, c,
    a + b, a + c, b + c, a + b + c
]

# Define edges of the unit cell
edges = [
    (0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4),
    (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)
]

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the unit cell edges
for edge in edges:
    p1, p2 = vertices[edge[0]], vertices[edge[1]]
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='black')

# (Optional) Plot atoms
# atom_positions = [np.array([0.5, 0.5, 0.5])]  # Example atom position
# ax.scatter([atom[0] for atom in atom_positions],
#            [atom[1] for atom in atom_positions],
#            [atom[2] for atom in atom_positions],
#            color='red', s=100)

# Set axis labels and limits
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.2, 1.2])
ax.set_zlim([-0.2, 1.2])

ax.set_box_aspect([1.0, 1.0, 1.0])

# Set title
ax.set_title('Unit Cell of Crystal')

# Show the plot
plt.show()