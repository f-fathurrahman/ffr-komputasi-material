import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. Define the potential and its derivatives
def morse_potential(r, De, a, r0):
    return De * ((1 - np.exp(-a*(r-r0)))**2 - 1)

def morse_force(r, De, a, r0):
    return -2*De*a*(np.exp(-a*(r-r0)) - np.exp(-2*a*(r-r0)))

def morse_force_constant(r, De, a, r0):
    return 2*De*a**2 * (2*np.exp(-2*a*(r-r0)) - np.exp(-a*(r-r0)))

# 2. Set realistic parameters (e.g., for a metal)
De = 2.0  # eV
a = 1.5   # Å⁻¹
r0 = 2.5  # Å (equilibrium distance)

# 3. Calculate force constants up to nth neighbor
max_neighbors = 3
force_constants = []
for n in range(1, max_neighbors+1):
    r = n*r0  # distance to nth neighbor
    phi = morse_force_constant(r, De, a, r0)
    force_constants.append(phi)
    print(f"Force constant for neighbor {n}: {phi:.4f} eV/Å²")

# 4. Build dynamical matrix
def dynamical_matrix(q, force_constants, r0, m):
    """Calculate dynamical matrix for a monoatomic chain"""
    D = 0j
    for n, phi in enumerate(force_constants, 1):
        r = n*r0
        D += 2*phi * (1 - np.cos(q*r)) / m
    return D

# 5. Calculate phonon dispersion
m = 26.98  # mass of aluminum in amu
q_points = np.linspace(-np.pi/r0, np.pi/r0, 100)
frequencies = []

for q in q_points:
    D = dynamical_matrix(q, force_constants, r0, m)
    # Convert from eV/Å²/amu to THz (1 eV/Å²/amu ≈ 15.6 THz²)
    frequencies.append(np.sqrt(np.abs(D)*15.6)/(2*np.pi))

# Plotting
plt.figure(figsize=(10,6))
plt.plot(q_points, frequencies, 'b-', linewidth=2)
plt.xlabel('Wave vector q (Å⁻¹)', fontsize=12)
plt.ylabel('Frequency (THz)', fontsize=12)
plt.title('Phonon Dispersion with Morse Potential', fontsize=14)
plt.grid(True)
plt.show()

