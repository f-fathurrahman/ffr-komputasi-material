import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ========================
# Parameters
# ========================
a = 5.43  # Silicon lattice constant (Å)
G_cutoff = 5  # Plane-wave cutoff (small for testing)
n_bands = 5  # Number of bands to compute
eta = 0.05  # Small broadening (Ry)
omega_max = 5.0  # Frequency cutoff (Ry)
n_omega = 20  # Frequency points (reduced for speed)

# High-symmetry k-point path (Γ-X-K-Γ)
k_points = np.array([
    [0.0, 0.0, 0.0],    # Γ
    [0.5, 0.0, 0.0],    # X
    [0.375, 0.375, 0.0], # K (typical for FCC)
    [0.0, 0.0, 0.0]     # Γ
])

# ========================
# Step 1: Generate G-vectors
# ========================
G_list = []
for h in range(-G_cutoff, G_cutoff+1):
    for k in range(-G_cutoff, G_cutoff+1):
        for l in range(-G_cutoff, G_cutoff+1):
            G = 2*np.pi/a * np.array([h,k,l])
            if np.linalg.norm(G) <= G_cutoff * 2*np.pi/a:
                G_list.append(G)
G_list = np.array(G_list)
N_G = len(G_list)
print("N_G = ", N_G)

# ========================
# Step 2: Model Potential (Empirical Pseudopotential)
# ========================
V_pseudo = np.zeros(N_G)
for i, G in enumerate(G_list):
    if np.linalg.norm(G) > 1e-8:  # Exclude G=0
        # Simple model: V(G) decreases with |G|
        V_pseudo[i] = -0.3 * np.exp(-0.5 * np.linalg.norm(G)**2)

# ========================
# Step 3: Solve Kohn-Sham Hamiltonian
# ========================
KS_bands = []
KS_wavefunctions = []
for k in k_points:
    H = np.zeros((N_G, N_G), dtype=complex)
    for i, Gi in enumerate(G_list):
        for j, Gj in enumerate(G_list):
            if i == j:
                H[i,j] = 0.5 * np.linalg.norm(k + Gi)**2  # Kinetic term
            # Potential term V(Gi-Gj)
            delta_G = Gi - Gj
            idx = np.argmin(np.linalg.norm(delta_G - G_list, axis=1))
            H[i,j] += V_pseudo[idx]
    
    eigvals, eigvecs = eigh(H)
    KS_bands.append(eigvals[:n_bands])  # Keep only lowest n_bands
    KS_wavefunctions.append(eigvecs[:, :n_bands])

KS_bands = np.array(KS_bands).T  # Shape: (n_bands, n_kpoints)
KS_wavefunctions = np.array(KS_wavefunctions)  # Shape: (n_kpoints, N_G, n_bands)

# ========================
# Step 4: GW Self-Energy Calculation
# ========================
def v_q_3D(q):
    """3D Coulomb interaction (regularized)"""
    return 4 * np.pi / (np.linalg.norm(q)**2 + 1e-10)

def compute_chi0(q, omega, ik, KS_bands, KS_wavefunctions, eta):
    """Compute χ₀(q,ω) at k-point ik"""
    chi0 = 0.0
    for n in range(n_bands):
        for m in range(n_bands):
            delta_e = KS_bands[m, ik] - KS_bands[n, ik]
            # Matrix element |⟨ψ_nk|e^{iq·r}|ψ_mk⟩|²
            psi_n = KS_wavefunctions[ik, :, n]
            psi_m = KS_wavefunctions[ik, :, m]
            M = np.abs(np.vdot(psi_n, np.exp(1j * np.dot(q, G_list.T)) * psi_m))**2
            chi0 += M * (1.0/(omega - delta_e + 1j*eta) - 1.0/(omega + delta_e + 1j*eta))
    return chi0

def W_q_3D(q, omega, ik, KS_bands, KS_wavefunctions, eta):
    """Compute W(q,ω) = v(q)/ε(q,ω)"""
    chi0 = compute_chi0(q, omega, ik, KS_bands, KS_wavefunctions, eta)
    return v_q_3D(q) / (1.0 - v_q_3D(q) * chi0)

def compute_self_energy_3D(ik, n, omega, KS_bands, KS_wavefunctions, eta):
    """Compute ⟨ψ_nk|Σ(ω)|ψ_nk⟩"""
    q = k_points[ik]
    Sigma = 0.0
    omega_grid = np.linspace(-omega_max, omega_max, n_omega)
    d_omega = omega_grid[1] - omega_grid[0]
    
    for omega_prime in omega_grid:
        # Green's function G(k, ω+ω')
        G = 0.0
        for m in range(n_bands):
            overlap = np.abs(np.vdot(KS_wavefunctions[ik, :, n], KS_wavefunctions[ik, :, m]))**2
            G += overlap / (omega + omega_prime - KS_bands[m, ik] + 1j*eta)
        
        # Screened interaction W(q, ω')
        W = W_q_3D(q, omega_prime, ik, KS_bands, KS_wavefunctions, eta)
        Sigma += G * W * np.exp(-1j*eta*omega_prime)
    
    Sigma *= 1j/(2*np.pi) * d_omega
    return np.real(Sigma)

# ========================
# Step 5: Apply GW Correction
# ========================
#GW_bands = np.zeros_like(KS_bands)
#for ik in range(len(k_points)):
#    for n in range(n_bands):
#        Sigma_nk = compute_self_energy_3D(ik, n, KS_bands[n, ik], KS_bands, KS_wavefunctions, eta)
#        GW_bands[n, ik] = KS_bands[n, ik] + Sigma_nk

# ========================
# Step 6: Plot Results
# ========================
plt.figure(figsize=(10, 6))
for band in range(n_bands):
    plt.plot(range(len(k_points)), KS_bands[band], 'b-', label='KS (DFT)' if band == 0 else "")
    plt.plot(range(len(k_points)), GW_bands[band], 'r--', label='GW' if band == 0 else "")
plt.xticks(range(len(k_points)), ["Γ", "X", "K", "Γ"])
plt.xlabel("k-path")
plt.ylabel("Energy (Ry)")
plt.title("3D Silicon Band Structure: KS vs GW")
plt.legend()
plt.grid(True)
plt.show()