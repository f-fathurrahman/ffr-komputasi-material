import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.fft import fft, ifft

# Constants
hbar = 1.0545718e-34  # Reduced Planck's constant (J·s)
e = 1.60217662e-19    # Electron charge (C)
kB = 1.38064852e-23   # Boltzmann constant (J/K)
m_e = 9.10938356e-31  # Electron mass (kg)

class Crystal1D:
    def __init__(self, a=5e-10, N=100, T=300, E_F=2.0, tau=1e-14):
        """
        Initialize 1D crystal model
        
        Parameters:
        a: lattice constant (m)
        N: number of unit cells
        T: temperature (K)
        E_F: Fermi energy (eV)
        tau: relaxation time (s)
        """
        self.a = a
        self.N = N
        self.T = T
        self.E_F = E_F * e  # Convert to Joules
        self.tau = tau
        
        # Generate reciprocal space (k-space)
        self.k = np.linspace(-np.pi/a, np.pi/a, N, endpoint=False)
        
        # Simple tight-binding dispersion relation E(k) = -2t*cos(ka)
        self.t = 1.6 * e  # Hopping parameter (J), ~1.6 eV
        self.Ek = -2 * self.t * np.cos(self.k * self.a)
        
        # Velocity v(k) = (1/hbar) * dE/dk
        self.vk = (2 * self.t * self.a / hbar) * np.sin(self.k * self.a)
        
    def fermi_dirac(self, E):
        """Fermi-Dirac distribution function"""
        return 1 / (np.exp((E - self.E_F) / (kB * self.T)) + 1)
    
    def dfdE(self, E):
        """Derivative of Fermi-Dirac distribution"""
        exp_term = np.exp((E - self.E_F) / (kB * self.T))
        return -exp_term / (kB * T * (exp_term + 1)**2)
    
    def compute_conductivity(self):
        """
        Compute electrical conductivity using Boltzmann Transport Equation
        
        Returns:
        sigma: electrical conductivity (S/m)
        """
        # Calculate -df/dE at each k-point
        df_dE = -self.dfdE(self.Ek)
        
        # Integrate using BTE in relaxation time approximation
        integrand = e**2 * self.tau * self.vk**2 * df_dE
        
        # 1D integration in k-space
        sigma = simps(integrand, self.k) / (np.pi * self.a)
        
        return sigma
    
    def plot_band_structure(self):
        """Plot the electronic band structure"""
        plt.figure(figsize=(10, 6))
        plt.plot(self.k, self.Ek / e, 'b-')  # Convert back to eV for plotting
        plt.axhline(self.E_F / e, color='r', linestyle='--', label='Fermi Energy')
        plt.xlabel('Wavevector k (1/m)', fontsize=12)
        plt.ylabel('Energy (eV)', fontsize=12)
        plt.title('1D Tight-Binding Band Structure', fontsize=14)
        plt.legend()
        plt.grid(True)
        plt.show()
    
    def plot_velocity(self):
        """Plot the electron velocity vs k"""
        plt.figure(figsize=(10, 6))
        plt.plot(self.k, self.vk, 'g-')
        plt.xlabel('Wavevector k (1/m)', fontsize=12)
        plt.ylabel('Electron velocity v(k) (m/s)', fontsize=12)
        plt.title('Electron Velocity in 1D Crystal', fontsize=14)
        plt.grid(True)
        plt.show()

# Example usage
if __name__ == "__main__":
    # Create a 1D crystal
    a = 5e-10  # Lattice constant (5 Å)
    N = 1000   # Number of k-points
    T = 300    # Temperature (K)
    E_F = 2.0  # Fermi energy (eV)
    tau = 1e-14  # Relaxation time (s)
    
    crystal = Crystal1D(a=a, N=N, T=T, E_F=E_F, tau=tau)
    
    # Plot band structure and velocity
    crystal.plot_band_structure()
    crystal.plot_velocity()
    
    # Compute conductivity
    sigma = crystal.compute_conductivity()
    print(f"Electrical conductivity: {sigma:.3f} S/m")
    
    # Temperature dependence study
    temps = np.linspace(50, 500, 50)
    conductivities = []
    
    for temp in temps:
        crystal.T = temp
        conductivities.append(crystal.compute_conductivity())
    
    plt.figure(figsize=(10, 6))
    plt.plot(temps, conductivities, 'r-')
    plt.xlabel('Temperature (K)', fontsize=12)
    plt.ylabel('Conductivity (S/m)', fontsize=12)
    plt.title('Temperature Dependence of Conductivity', fontsize=14)
    plt.grid(True)
    plt.show()


