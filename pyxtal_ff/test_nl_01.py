import numpy as np
from itertools import combinations

def calculate_neighbor_list(positions, cutoff, box, pbc=True):
    """
    Calculate neighbor list for a set of atoms with periodic boundary conditions.
    
    Args:
        positions (np.array): Nx3 array of atomic positions
        cutoff (float): Distance cutoff for neighbor consideration
        box (np.array): 3-element array of box dimensions [Lx, Ly, Lz]
        pbc (bool or list): Periodic boundary conditions flag (True/False or [x,y,z])
    
    Returns:
        list: List of tuples (i, j, distance) for each neighbor pair
    """
    if isinstance(pbc, bool):
        pbc = [pbc, pbc, pbc]
    
    N = len(positions)
    neighbor_list = []
    
    # Check all unique pairs
    for i, j in combinations(range(N), 2):
        rij = positions[j] - positions[i]
        
        # Apply minimum image convention for periodic directions
        for dim in range(3):
            if pbc[dim] and box[dim] > 0:
                rij[dim] -= round(rij[dim] / box[dim]) * box[dim]
        
        distance = np.linalg.norm(rij)
        
        if distance <= cutoff:
            neighbor_list.append((i, j, distance))
    
    return neighbor_list

def calculate_neighbor_list_cell(positions, cutoff, box, pbc=True):
    """
    More efficient version using cell lists for large systems.
    """
    if isinstance(pbc, bool):
        pbc = [pbc, pbc, pbc]
    
    N = len(positions)
    neighbor_list = []
    
    # Determine optimal number of cells
    min_box_dim = min(box[pbc]) if any(pbc) else cutoff
    n_cells = max(3, int(min_box_dim / cutoff))
    
    # Create cell list
    cell_size = box / n_cells
    cell_list = [[] for _ in range(n_cells**3)]
    
    # Assign atoms to cells
    for i, pos in enumerate(positions):
        cell_idx = tuple((pos / cell_size).astype(int) % n_cells)
        linear_idx = np.ravel_multi_index(cell_idx, (n_cells, n_cells, n_cells))
        cell_list[linear_idx].append(i)
    
    # Check neighboring cells
    for cell_linear in range(len(cell_list)):
        cell_idx = np.unravel_index(cell_linear, (n_cells, n_cells, n_cells))
        
        # Get all atoms in current cell
        atoms_in_cell = cell_list[cell_linear]
        
        # Get neighboring cells (including self)
        for dx, dy, dz in [(0,0,0), (1,0,0), (-1,0,0), 
                           (0,1,0), (0,-1,0), (1,1,0), 
                           (-1,-1,0), (1,-1,0), (-1,1,0),
                           # Add more offsets for 3D if needed
                           ]:
            neighbor_idx = ((cell_idx[0] + dx) % n_cells,
                           (cell_idx[1] + dy) % n_cells,
                           (cell_idx[2] + dz) % n_cells)
            neighbor_linear = np.ravel_multi_index(neighbor_idx, (n_cells, n_cells, n_cells))
            
            for i in atoms_in_cell:
                for j in cell_list[neighbor_linear]:
                    if i >= j:  # Avoid duplicate pairs
                        continue
                        
                    rij = positions[j] - positions[i]
                    
                    # Apply minimum image convention
                    for dim in range(3):
                        if pbc[dim] and box[dim] > 0:
                            rij[dim] -= round(rij[dim] / box[dim]) * box[dim]
                    
                    distance = np.linalg.norm(rij)
                    
                    if distance <= cutoff:
                        neighbor_list.append((i, j, distance))
    
    return neighbor_list


positions = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [9.0, 0.0, 0.0]  # This should be a neighbor of [0,0,0] with PBC
])

box = np.array([10.0, 10.0, 10.0])
cutoff = 1.5

# Simple version (good for small systems)
neighbors = calculate_neighbor_list(positions, cutoff, box)
print("Simple neighbor list:")
for i, j, d in neighbors:
    print(f"Atoms {i}-{j}: distance {d:.3f}")

# Cell list version (better for large systems)
neighbors_cell = calculate_neighbor_list_cell(positions, cutoff, box)
print("\nCell list neighbor list:")
for i, j, d in neighbors_cell:
    print(f"Atoms {i}-{j}: distance {d:.3f}")


