
import numpy as np
from typing import Literal, Optional


from dataclasses import dataclass
from .occupation_states import OccupationInfo
from ..scf.density import DensityData
from ..scf.energy import EnergyComponents
from ..scf.driver import IntermediateInfo
from ..mesh.operators import GridData


ResultType = Literal["SCF", "forward"]


@dataclass
class AtomicDFTSolverResult:
    """
    Result from AtomicDFTSolver
    """
    # Result type, SCF or forward
    result_type               : ResultType                 # SCF or forward pass

    # electron density related information
    rho                       : np.ndarray                 # Electron density
    rho_nlcc                  : np.ndarray                 # Non-linear core correction density
    density_data              : DensityData                # Density data

    # energy related information
    energy                    : float                      # Total energy
    energy_components         : EnergyComponents           # Energy components
    

    # spectrum related information
    eigen_energies            : Optional[np.ndarray]       # Kohn-Sham eigenvalues
    orbitals                  : np.ndarray                 # Kohn-Sham orbitals
    orbitals_on_uniform_grid  : np.ndarray                 # Kohn-Sham orbitals on uniform grid
    full_eigen_energies       : Optional[np.ndarray]       # Full eigenvalues
    full_orbitals             : Optional[np.ndarray]       # Full orbitals
    full_l_terms              : Optional[np.ndarray]       # Full l terms
    occupation_info           : OccupationInfo             # Occupation info

    # grid related information
    quadrature_nodes          : np.ndarray                 # Quadrature nodes
    quadrature_weights        : np.ndarray                 # Quadrature weights
    uniform_grid              : np.ndarray                 # Uniform grid
    grid_data                 : GridData                   # Grid data

    # Potential and energy density related information
    v_x_local                 : np.ndarray                 # Local XC potential   
    v_c_local                 : np.ndarray                 # Local XC potential
    v_x_local_on_uniform_grid : np.ndarray                 # Local XC potential on uniform grid
    v_c_local_on_uniform_grid : np.ndarray                 # Local XC potential on uniform grid
    e_x_local                 : Optional[np.ndarray]       # Local XC energy density
    e_c_local                 : Optional[np.ndarray]       # Local XC energy density
    e_x_local_on_uniform_grid : Optional[np.ndarray]       # Local XC energy density on uniform grid
    e_c_local_on_uniform_grid : Optional[np.ndarray]       # Local XC energy density on uniform grid

    # SCF related information
    converged                 : Optional[bool]             # Whether the SCF converged
    iterations                : Optional[int]              # Number of SCF iterations
    rho_residual              : Optional[float]            # Density residual
    intermediate_info         : Optional[IntermediateInfo] # Intermediate information from SCF iterations


    def __post_init__(self):
        pass


    
