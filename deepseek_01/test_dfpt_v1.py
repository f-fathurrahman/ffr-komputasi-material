import numpy as np
from scipy.fft import fft, ifft
from scipy.linalg import eigh

# Parameters
a = 5.0  # lattice constant
Nk = 32  # number of k-points
Ng = 32  # number of plane waves
n_atoms = 2  # atoms per unit cell
atom_positions = np.array([0.4, 0.6]) * a  # positions within unit cell
Z = [1, 1]  # atomic numbers (simplified)

# Reciprocal space grid
G = np.fft.fftfreq(Ng, d=a/Ng) * 2*np.pi
k_points = np.linspace(0, 2*np.pi/a, Nk, endpoint=False)

# Pseudopotential (simplified)
#def v_ps(G, Z):
#    return -Z * np.exp(-(G*0.2)**2) / (G**2 + 1e-6)

# Reciprocal space form (for plane wave calculations)
def v_ps(G, Z, rc=1.5):
    """Fourier transform of the softened potential"""
    return -Z * np.exp(-rc*abs(G)) * (1 + abs(G)/2)

# Construct Hamiltonian for each k-point
def solve_ks(k, n_bands=4):
    H = np.zeros((Ng, Ng), dtype=complex)
    
    # Kinetic energy
    for i, gi in enumerate(G):
        H[i,i] = 0.5 * (k + gi)**2
        
    # Local potential
    for i, gi in enumerate(G):
        for j, gj in enumerate(G):
            G_diff = gi - gj
            for pos, Z_i in zip(atom_positions, Z):
                H[i,j] += v_ps(G_diff, Z_i) * np.exp(-1j*G_diff*pos) / a
    
    # Solve
    eps, psi = eigh(H)
    return eps[:n_bands], psi[:, :n_bands]



def calculate_phonons():
    # First, calculate ground state
    n_bands = 4
    n_k = len(k_points)
    psi_k = np.zeros((n_k, Ng, n_bands), dtype=complex)
    eps_k = np.zeros((n_k, n_bands))
    
    for ik, k in enumerate(k_points):
        eps_k[ik], psi_k[ik] = solve_ks(k, n_bands)
    

    # Now DFPT for phonons
    dynamical_matrix = np.zeros((n_atoms, n_atoms))
    
    for alpha in range(n_atoms):  # displaced atom
        for beta in range(n_atoms):  # force on atom
            D_alpha_beta = 0.0
            
            # Calculate the derivative of the potential
            for ik, k in enumerate(k_points):
                for n in range(n_bands):  # occupied bands
                    # Perturbation potential derivative
                    dV_psi = np.zeros(Ng, dtype=complex)
                    for i, gi in enumerate(G):
                        for j, gj in enumerate(G):
                            G_diff = gi - gj
                            dV = -1j * G_diff * v_ps(G_diff, Z[beta]) * np.exp(-1j*G_diff*atom_positions[beta])
                            dV_psi[i] += dV * psi_k[ik, j, n] / a
                    
                    # Modified response function calculation
                    response = np.zeros(Ng, dtype=complex)
                    for m in range(n_bands):
                        if m == n or abs(eps_k[ik, n] - eps_k[ik, m]) < 1e-10:
                            continue  # Skip degenerate or same states
                        numerator = np.vdot(psi_k[ik, :, m], dV_psi)
                        denominator = eps_k[ik, n] - eps_k[ik, m]
                        response += numerator / denominator * psi_k[ik, :, m]


                    # Add to dynamical matrix
                    for i, gi in enumerate(G):
                        for j, gj in enumerate(G):
                            G_diff = gi - gj
                            term = v_ps(G_diff, Z[alpha]) * np.exp(-1j*G_diff*atom_positions[alpha])
                            term *= np.conj(psi_k[ik, i, n]) * response[j]
                            D_alpha_beta += -2 * np.real(term) / n_k  # factor 2 for spin
            
            dynamical_matrix[alpha, beta] = D_alpha_beta
    
    # Add the ionic contribution (Ewald sum in 1D)
    for alpha in range(n_atoms):
        for beta in range(n_atoms):
            if alpha == beta:
                for gamma in range(n_atoms):
                    if gamma != alpha:
                        r = atom_positions[alpha] - atom_positions[gamma]
                        dynamical_matrix[alpha, alpha] += Z[alpha]*Z[gamma] * 2*np.pi/a * (1/a - np.abs(r)/a**2)
            else:
                r = atom_positions[alpha] - atom_positions[beta]
                dynamical_matrix[alpha, beta] -= Z[alpha]*Z[beta] * 2*np.pi/a * (1/a - np.abs(r)/a**2)
    
    return dynamical_matrix



def phonon_band_structure():
    n_q = 20
    q_points = np.linspace(0, 2*np.pi/a, n_q)
    frequencies = np.zeros((n_q, n_atoms))
    
    for iq, q in enumerate(q_points):
        # Calculate dynamical matrix at q
        D_q = np.zeros((n_atoms, n_atoms), dtype=complex)
        
        # This would involve modifying the DFPT calculation to include q
        # For simplicity, we'll just use the Gamma point result here
        D_q = calculate_phonons()
        
        # Diagonalize to get frequencies
        eigvals = np.linalg.eigvalsh(D_q)
        frequencies[iq] = np.sqrt(np.abs(eigvals)) * np.sign(eigvals)

        print(f"DONE calc for iq = {iq} q = {q}")
    
    return q_points, frequencies




import matplotlib.pyplot as plt

q, freq = phonon_band_structure()

plt.figure(figsize=(8,6))
for i in range(n_atoms):
    plt.plot(q, freq[:, i])

plt.xlabel('q (2π/a)')
plt.ylabel('Frequency (arb. units)')
plt.title('1D Phonon Band Structure')
plt.show()



