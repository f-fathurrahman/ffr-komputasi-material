import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ========================
# Parameters
# ========================
a = 5.43  # Silicon lattice constant (Å)
G_cutoff = 2  # Reduced for testing
n_bands = 4  # Number of bands to compute
n_kpoints_per_segment = 10  # k-points between high-symmetry points

# ========================
# FCC Lattice Setup (corrected)
# ========================
a1 = np.array([0.0, 0.5, 0.5]) * a
a2 = np.array([0.5, 0.0, 0.5]) * a
a3 = np.array([0.5, 0.5, 0.0]) * a

# Reciprocal lattice vectors
b1 = 2*np.pi*np.cross(a2, a3)/np.dot(a1, np.cross(a2, a3))
b2 = 2*np.pi*np.cross(a3, a1)/np.dot(a1, np.cross(a2, a3))
b3 = 2*np.pi*np.cross(a1, a2)/np.dot(a1, np.cross(a2, a3))

# Generate G-vectors
max_n = int(np.ceil(G_cutoff))
G_list = []
for h in range(-max_n, max_n+1):
    for k in range(-max_n, max_n+1):
        for l in range(-max_n, max_n+1):
            G = h*b1 + k*b2 + l*b3
            if np.linalg.norm(G) <= G_cutoff*2*np.pi/a:
                G_list.append(G)
G_list = np.array(G_list)
N_G = len(G_list)

# Structure factor for diamond (2-atom basis)
def structure_factor(G):
    tau1 = np.array([0.0, 0.0, 0.0])
    tau2 = np.array([0.25, 0.25, 0.25]) * a
    return 1 + np.exp(-1j*np.dot(G, tau2))

# Empirical pseudopotential
def V_pseudo(G):
    G_norm = np.linalg.norm(G)
    if G_norm < 1e-8:
        return 0.0
    return -3.0 * np.exp(-0.5*(G_norm*a/(2*np.pi))**2) * np.real(structure_factor(G))

# ========================
# k-Path Generation (corrected)
# ========================
high_sym_kpoints = {
    'Γ': np.array([0.0, 0.0, 0.0]),
    'X': np.array([0.5, 0.0, 0.5]),
    'L': np.array([0.5, 0.5, 0.5]),
    'Γ2': np.array([0.0, 0.0, 0.0])
}

k_path_segments = ['Γ', 'X', 'L', 'Γ2']

k_points = []
label_positions = [0]
for i in range(len(k_path_segments)-1):
    start = high_sym_kpoints[k_path_segments[i]]
    end = high_sym_kpoints[k_path_segments[i+1]]
    for j in range(n_kpoints_per_segment):
        alpha = j/(n_kpoints_per_segment-1)
        k_points.append((1-alpha)*start + alpha*end)
    label_positions.append(len(k_points)-1)
k_points = np.array(k_points)

# ========================
# Kohn-Sham Solution (corrected array shapes)
# ========================
KS_bands = np.zeros((len(k_points), n_bands))
KS_wavefunctions = np.zeros((len(k_points), N_G, n_bands), dtype=complex)

for ik, k in enumerate(k_points):
    H = np.zeros((N_G, N_G), dtype=complex)
    for i, Gi in enumerate(G_list):
        for j, Gj in enumerate(G_list):
            if i == j:
                H[i,j] = 0.5 * np.linalg.norm(k + Gi)**2
            H[i,j] += V_pseudo(Gi - Gj)
    
    eigvals, eigvecs = eigh(H)
    KS_bands[ik, :] = eigvals[:n_bands]
    KS_wavefunctions[ik, :, :] = eigvecs[:, :n_bands]

# Transpose for plotting: shape (n_bands, n_kpoints)
KS_bands = KS_bands.T

# ========================
# Plotting (fixed)
# ========================
plt.figure(figsize=(12, 6))
x_axis = np.arange(len(k_points))

for band in range(n_bands):
    plt.plot(x_axis, KS_bands[band], 'b-', lw=1.5, alpha=0.7, 
             label=f'Band {band+1}')

plt.xticks(label_positions, k_path_segments)
for pos in label_positions:
    plt.axvline(pos, color='k', linestyle='--', alpha=0.2)

plt.xlabel("k-path")
plt.ylabel("Energy (Ry)")
plt.title("Silicon Band Structure (DFT)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()