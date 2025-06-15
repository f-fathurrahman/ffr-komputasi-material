import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# Parameters
a = 1.0                # Lattice constant
V0 = 0.5               # Potential strength
G_cutoff = 25          # Reduced for faster testing
k_points = 100         # Reduced k-points
n_bands = 4            # Bands to consider
eta = 0.05             # Broadening
omega_max = 5.0        # Frequency cutoff
n_omega = 50           # Frequency points

# Reciprocal lattice vectors
G_vec = np.array([n * (2 * np.pi / a) for n in range(-G_cutoff, G_cutoff + 1)])
N_G = len(G_vec)
print("N_G = ", N_G)


# Brillouin zone sampling
k_vec = np.linspace(-np.pi/a, np.pi/a, k_points)

# --- Step 1: Solve KS Hamiltonian ---
def solve_KS_bands(k_vec, G_vec, V0, a):
    bands = []
    wavefunctions = []
    for k in k_vec:
        H = np.zeros((N_G, N_G), dtype=complex)
        for i, G in enumerate(G_vec):
            for j, G_prime in enumerate(G_vec):
                if i == j:
                    H[i,j] = 0.5 * (k + G)**2
                if abs(G - G_prime) == 2 * np.pi / a:
                    H[i,j] += V0 / 2
        eigvals, eigvecs = eigh(H)
        bands.append(eigvals[:n_bands])  # Keep only lowest n_bands
        wavefunctions.append(eigvecs[:, :n_bands])
    return np.array(bands).T, np.array(wavefunctions)

KS_bands, KS_wavefunctions = solve_KS_bands(k_vec, G_vec, V0, a)

# --- Step 2: Helper functions for χ₀, W ---
def v_q(q):
    """1D Coulomb interaction (regularized)"""
    return 2 * np.log(np.abs(q) + 1e-10)  # Prevents divergence

def compute_chi0(q, omega, ik, KS_bands, KS_wavefunctions, eta):
    """χ₀(q,ω) at k-point ik"""
    chi0 = 0.0
    for n in range(n_bands):
        for m in range(n_bands):
            delta_e = KS_bands[m, ik] - KS_bands[n, ik]
            M = np.abs(np.vdot(KS_wavefunctions[ik, :, n], KS_wavefunctions[ik, :, m]))**2
            chi0 += M * (1.0/(omega - delta_e + 1j*eta) - 1.0/(omega + delta_e + 1j*eta))
    return chi0

def W_q(q, omega, ik, KS_bands, KS_wavefunctions, eta):
    """W(q,ω) = v(q) / ε(q,ω)"""
    chi0 = compute_chi0(q, omega, ik, KS_bands, KS_wavefunctions, eta)
    return v_q(q) / (1.0 - v_q(q) * chi0)

# --- Step 3: Compute Σ_nk ---
def compute_self_energy(ik, n, omega, KS_bands, KS_wavefunctions, eta):
    """⟨ψ_nk|Σ(ω)|ψ_nk⟩"""
    q = k_vec[ik]
    Sigma = 0.0
    omega_grid = np.linspace(-omega_max, omega_max, n_omega)
    d_omega = omega_grid[1] - omega_grid[0]
    
    for omega_prime in omega_grid:
        G = 0.0
        for m in range(n_bands):
            overlap = np.abs(np.vdot(KS_wavefunctions[ik, :, n], KS_wavefunctions[ik, :, m]))**2
            G += overlap / (omega + omega_prime - KS_bands[m, ik] + 1j*eta)
        W = W_q(q, omega_prime, ik, KS_bands, KS_wavefunctions, eta)
        Sigma += G * W * np.exp(-1j*eta*omega_prime)
    
    Sigma *= 1j/(2*np.pi) * d_omega
    return np.real(Sigma)

# --- Step 4: GW Correction ---
GW_bands = np.zeros_like(KS_bands)
for ik in range(len(k_vec)):
    for n in range(n_bands):
        Sigma_nk = compute_self_energy(ik, n, KS_bands[n, ik], KS_bands, KS_wavefunctions, eta)
        GW_bands[n, ik] = KS_bands[n, ik] + Sigma_nk

# --- Plotting ---
plt.figure(figsize=(8, 5))
for band in range(n_bands):
    plt.plot(k_vec, KS_bands[band], 'b-', label='KS (DFT)' if band == 0 else "")
    plt.plot(k_vec, GW_bands[band], 'r--', label='GW' if band == 0 else "")
plt.xlabel('k (1/a)')
plt.ylabel('Energy (eV)')
plt.title('1D Band Structure: KS vs GW (Corrected)')
plt.legend()
plt.grid(True)
plt.show()