import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

class PhononCalculator1D:
    def __init__(self, basis, masses, De, a, r0, max_neighbors=3):
        """
        Generalized 1D phonon calculator
        
        Parameters:
        - basis: List of atomic positions within one unit cell (in fractional coordinates 0-1)
        - masses: List of atomic masses (in amu)
        - De: Morse potential depth (eV)
        - a: Morse potential width parameter (Å⁻¹)
        - r0: Equilibrium bond length (Å)
        - max_neighbors: Number of neighbor shells to consider
        """
        self.basis = np.array(basis)
        self.masses = np.array(masses)
        self.De = De
        self.a = a
        self.r0 = r0
        self.max_neighbors = max_neighbors
        self.natoms = len(basis)
        
        # Set lattice constant - now handles monoatomic case
        self.a_lat = r0  # Default to r0 for monoatomic
        if self.natoms > 1:
            self.a_lat = self._calculate_lattice_constant()
    
    def _calculate_lattice_constant(self):
        """Determine the smallest lattice constant that tiles all basis positions"""
        sorted_basis = np.sort(self.basis)
        diffs = np.diff(sorted_basis)
        diffs = np.append(diffs, 1 + sorted_basis[0] - sorted_basis[-1])  # Wrap around
        return np.min(diffs[diffs > 0]) * self.r0 / np.min(diffs[diffs > 0])
    
    def morse_force_constant(self, r):
        """Calculate force constant from Morse potential at distance r"""
        return 2*self.De*self.a**2 * (2*np.exp(-2*self.a*(r-self.r0)) - np.exp(-self.a*(r-self.r0)))
    
    def build_force_constants(self):
        """Build force constant matrix considering all atom pairs within cutoff"""
        self.phi = np.zeros((self.natoms, self.natoms, self.max_neighbors))
        
        for i in range(self.natoms):
            for j in range(self.natoms):
                for n in range(self.max_neighbors):
                    # Calculate distance between atom i in home cell and atom j in nth neighbor cell
                    if self.natoms == 1:
                        # Monoatomic case - simple chain
                        r_ij = (n+1)*self.r0
                    else:
                        # Polyatomic case
                        r_ij = np.abs(self.basis[j] - self.basis[i] + (n+1)) * self.a_lat
                    self.phi[i,j,n] = self.morse_force_constant(r_ij)
        
        # Apply acoustic sum rule
        self._enforce_sum_rules()
    
    def _enforce_sum_rules(self):
        """Enforce translational invariance (sum of force constants = 0)"""
        for i in range(self.natoms):
            for n in range(self.max_neighbors):
                total = np.sum(self.phi[i,:,n])
                self.phi[i,i,n] -= total
    
    def dynamical_matrix(self, q):
        """Build dynamical matrix for wavevector q"""
        D = np.zeros((self.natoms, self.natoms), dtype=complex)
        
        for i in range(self.natoms):
            for j in range(self.natoms):
                for n in range(self.max_neighbors):
                    # Phase factor for (n+1)th neighbor
                    phase = np.exp(-1j*q*(n+1)*self.a_lat)
                    # Mass-weighted force constant
                    D[i,j] += self.phi[i,j,n] * phase / np.sqrt(self.masses[i]*self.masses[j])
        
        # Make Hermitian by averaging with conjugate transpose
        D = 0.5*(D + D.conj().T)
        return D
    
    def calculate_phonons(self, q_points):
        """Calculate phonon frequencies across q-points"""
        frequencies = []
        
        for q in q_points:
            D = self.dynamical_matrix(q)
            eigvals = eigh(D, eigvals_only=True)  # Returns sorted eigenvalues
            # Convert to THz: 1 eV/Å²/amu ≈ 15.6 THz²
            frequencies.append(np.sqrt(np.abs(eigvals)*15.6)/(2*np.pi))
        
        return np.array(frequencies)
    
    def plot_dispersion(self, q_points, frequencies, ylim=None):
        """Plot phonon dispersion relations"""
        plt.figure(figsize=(10,6))
        
        for band in range(frequencies.shape[1]):
            plt.plot(q_points, frequencies[:,band], 'b-')
        
        plt.xlabel('Wave vector q (Å⁻¹)', fontsize=12)
        plt.ylabel('Frequency (THz)', fontsize=12)
        plt.title(f'Phonon Dispersion (N={self.natoms} atoms/unit cell)', fontsize=14)
        if ylim is not None:
            plt.ylim(ylim)
        plt.grid(True)
        plt.show()

# Test Cases

# 1. Monoatomic chain (now works)
print("Monoatomic chain:")
mono_calc = PhononCalculator1D(
    basis=[0],          # One atom at position 0
    masses=[26.98],     # Aluminum mass
    De=2.0, a=1.5, r0=2.5
)
mono_calc.build_force_constants()
q_points = np.linspace(-np.pi/mono_calc.r0, np.pi/mono_calc.r0, 100)
freqs = mono_calc.calculate_phonons(q_points)
mono_calc.plot_dispersion(q_points, freqs, ylim=[0, 10])

# 2. Diatomic chain
print("\nDiatomic chain:")
di_calc = PhononCalculator1D(
    basis=[0, 0.5],     # Two atoms at 0 and 0.5
    masses=[26.98, 92.91],  # Al and Nb masses
    De=2.0, a=1.5, r0=2.5
)
di_calc.build_force_constants()
q_points = np.linspace(-np.pi/di_calc.a_lat, np.pi/di_calc.a_lat, 100)
freqs = di_calc.calculate_phonons(q_points)
di_calc.plot_dispersion(q_points, freqs, ylim=[0, 10])

# 3. Three-atom chain
print("\nThree-atom chain:")
tri_calc = PhononCalculator1D(
    basis=[0, 0.3, 0.7],
    masses=[12.0, 16.0, 12.0],
    De=3.0, a=1.8, r0=2.0
)
tri_calc.build_force_constants()
q_points = np.linspace(-np.pi/tri_calc.a_lat, np.pi/tri_calc.a_lat, 100)
freqs = tri_calc.calculate_phonons(q_points)
tri_calc.plot_dispersion(q_points, freqs, ylim=[0, 15])