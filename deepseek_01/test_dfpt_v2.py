import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# Parameters
a = 5.0  # Lattice constant
n_atoms = 2  # Diatomic chain
atom_positions = np.array([0.2, 0.7]) * a
masses = [1.0, 1.5]  # Different masses for dispersion
Z = [1, 2]  # Different atomic numbers

# Reciprocal space setup
Nk = 50  # k-points for electronic bands
Ng = 32  # Plane waves
G = np.fft.fftfreq(Ng, d=a/Ng) * 2*np.pi
k_path = np.linspace(0, 2*np.pi/a, Nk, endpoint=False)

# Potential functions
def v_ps_g(G, Z, rc=0.3):
    return -Z * np.exp(-rc*np.abs(G)) * (1 + rc*np.abs(G))/2

def build_hamiltonian(k):
    H = np.zeros((Ng, Ng), dtype=complex)
    
    # Kinetic energy
    for i, gi in enumerate(G):
        H[i,i] = 0.5 * (k + gi)**2
    
    # Potential energy
    for i, gi in enumerate(G):
        for j, gj in enumerate(G):
            G_diff = gi - gj
            for pos, Z_i in zip(atom_positions, Z):
                H[i,j] += v_ps_g(G_diff, Z_i) * np.exp(-1j*G_diff*pos) / a
    return H

def calculate_electronic_bands(n_bands=4):
    eigenvalues = np.zeros((len(k_path), n_bands))
    
    for ik, k in enumerate(k_path):
        eps, _ = eigh(build_hamiltonian(k))
        eigenvalues[ik] = eps[:n_bands]  # Only store lowest n_bands
    
    return eigenvalues

def calculate_phonon_bands():
    n_q = 30
    q_path = np.linspace(0, np.pi/a, n_q)
    frequencies = np.zeros((n_q, n_atoms))
    n_bands = 4  # Number of bands to consider
    
    # Store wavefunctions for occupied states
    psi_k = np.zeros((len(k_path), Ng, n_bands), dtype=complex)
    eps_k = np.zeros((len(k_path), n_bands))
    
    for ik, k in enumerate(k_path):
        eps_all, psi_all = eigh(build_hamiltonian(k))
        eps_k[ik] = eps_all[:n_bands]  # Only store lowest n_bands
        psi_k[ik] = psi_all[:, :n_bands]  # Only store lowest n_bands
    
    for iq, q in enumerate(q_path):
        D = np.zeros((n_atoms, n_atoms), dtype=complex)
        
        for alpha in range(n_atoms):
            for beta in range(n_atoms):
                for ik, k in enumerate(k_path):
                    for n in range(n_bands//2):  # Only occupied states
                        # Build perturbation potential
                        dV_beta = np.zeros(Ng, dtype=complex)
                        for i, gi in enumerate(G):
                            for j, gj in enumerate(G):
                                G_diff = gi - gj
                                dV = -1j*G_diff * v_ps_g(G_diff, Z[beta]) * np.exp(-1j*G_diff*atom_positions[beta])
                                dV_beta[i] += dV * psi_k[ik,j,n] / a
                        
                        # Calculate response
                        response = np.zeros(Ng, dtype=complex)
                        for m in range(n_bands):
                            if m == n: continue
                            numerator = np.vdot(psi_k[ik,:,m], dV_beta)
                            denominator = eps_k[ik,n] - eps_k[ik,m]
                            if abs(denominator) > 1e-6:
                                response += numerator / denominator * psi_k[ik,:,m]
                        
                        # Add to dynamical matrix
                        phase = np.exp(1j*q*(atom_positions[alpha]-atom_positions[beta]))
                        for i, gi in enumerate(G):
                            for j, gj in enumerate(G):
                                G_diff = gi - gj
                                term = -1j*G_diff * v_ps_g(G_diff, Z[alpha]) * np.exp(-1j*G_diff*atom_positions[alpha])
                                D[alpha,beta] += 2 * phase * np.conj(psi_k[ik,i,n]) * response[j] / len(k_path)
        
        # Add short-range repulsion
        for alpha in range(n_atoms):
            for beta in range(n_atoms):
                r = (atom_positions[alpha] - atom_positions[beta]) % a
                r = min(r, a-r)
                D[alpha,beta] += 10.0 * np.exp(-(r/0.8)**2) * (1 - 2*(r/0.8)**2) / np.sqrt(masses[alpha]*masses[beta])
        
        # Diagonalize dynamical matrix
        eigvals = np.linalg.eigvalsh(D)
        frequencies[iq] = np.sign(eigvals) * np.sqrt(np.abs(eigvals))
    
    return q_path, frequencies

def plot_bands():
    plt.figure(figsize=(12, 8))
    
    # Electronic bands
    plt.subplot(2, 1, 1)
    e_bands = calculate_electronic_bands()
    for band in range(e_bands.shape[1]):
        plt.plot(k_path*a/np.pi, e_bands[:, band], 'b-')
    plt.ylabel('Energy (eV)')
    plt.title('Electronic Bands')
    plt.xticks([0, 1, 2], ['Γ', 'X', '2π/a'])
    plt.grid(True)
    
    # Phonon bands
    plt.subplot(2, 1, 2)
    q_points, phonon_bands = calculate_phonon_bands()
    for band in range(phonon_bands.shape[1]):
        plt.plot(q_points*a/np.pi, phonon_bands[:, band], 'r-')
    plt.xlabel('Wavevector (π/a)')
    plt.ylabel('Frequency (THz)')
    plt.title('Phonon Dispersion')
    plt.xticks([0, 0.5, 1], ['Γ', 'X', 'π/a'])
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

plot_bands()