"""
Response function calculator for Kohn-Sham DFT
Computes response function (chi_0) and dielectric matrix for preconditioning
"""

from __future__ import annotations
import numpy as np

from ..mesh.operators import GridData, RadialOperatorsBuilder
from ..utils.occupation_states import OccupationInfo

# Error messages
GRID_DATA_TYPE_ERROR_MESSAGE = \
    "parameter 'grid_data' must be an instance of GridData, get type {} instead."
OCCUPATION_INFO_TYPE_ERROR_MESSAGE = \
    "parameter 'occupation_info' must be an instance of OccupationInfo, get type {} instead."
OPS_BUILDER_TYPE_ERROR_MESSAGE = \
    "parameter 'ops_builder' must be an instance of RadialOperatorsBuilder, get type {} instead."
FULL_EIGENVALUES_TYPE_ERROR_MESSAGE = \
    "parameter 'full_eigenvalues' must be a numpy array, get type {} instead."
FULL_ORBITALS_TYPE_ERROR_MESSAGE = \
    "parameter 'full_orbitals' must be a numpy array, get type {} instead."
FULL_L_TERMS_TYPE_ERROR_MESSAGE = \
    "parameter 'full_l_terms' must be a numpy array, get type {} instead."


class ResponseCalculator:
    """
    Computes response function (chi_0) and dielectric matrix from Kohn-Sham eigenstates.
    
    The response function describes how the electron density responds to changes in the
    external potential. The dielectric matrix is used as a preconditioner for Pulay mixing.
    
    For spherically symmetric atoms, the response function is computed using:
        χ₀(r, r') = Σ_{i,occ} Σ_{a,unocc} f_i (2l_i + 1) 
                    × [ψ_i(r) G_{ia}(r, r') ψ_i(r')] / (ε_i - ε_a)
    
    where G_{ia} is the Green's function block for transitions from occupied state i
    to unoccupied state a.
    """
    
    def __init__(
        self,
        occupation_info : OccupationInfo,
        ops_builder     : RadialOperatorsBuilder,
    ):
        """
        Parameters
        ----------
        occupation_info : OccupationInfo
            Occupation numbers and quantum numbers
        ops_builder : RadialOperatorsBuilder
            Operators builder for interpolation and differentiation
        """
        
        grid_data = ops_builder.grid_data

        assert isinstance(grid_data, GridData), \
            GRID_DATA_TYPE_ERROR_MESSAGE.format(type(grid_data))
        assert isinstance(occupation_info, OccupationInfo), \
            OCCUPATION_INFO_TYPE_ERROR_MESSAGE.format(type(occupation_info))
        assert isinstance(ops_builder, RadialOperatorsBuilder), \
            OPS_BUILDER_TYPE_ERROR_MESSAGE.format(type(ops_builder))
        
        self.quadrature_nodes   = grid_data.quadrature_nodes
        self.quadrature_weights = grid_data.quadrature_weights
        self.occupation_info    = occupation_info
        self.ops_builder        = ops_builder
        
        # Extract occupation information
        self.occ_l_values = occupation_info.l_values
        self.occ_n_values = occupation_info.n_values
        self.occupations  = occupation_info.occupations
        self.n_quad       = len(self.quadrature_nodes)
        self.n_interior   = len(ops_builder.physical_nodes) - 2
    

    def compute_chi_0_kernel(
        self,
        full_eigenvalues : np.ndarray,
        full_orbitals    : np.ndarray,
        full_l_terms     : np.ndarray,
    ) -> np.ndarray:
        """
        Compute the response function kernel χ₀(r, r').
        Here, we assume that the input full_orbitals is on the quadrature grid.

        Parameters
        ----------
        full_eigenvalues : np.ndarray
            Full eigenvalues (occupied + unoccupied), shape (n_states,)
        full_orbitals : np.ndarray
            Full orbitals at quadrature nodes, shape (n_states, n_quad)
        full_l_terms : np.ndarray
            Angular momentum quantum numbers, shape (n_states,)
        
        Returns
        -------
        chi_0_kernel : np.ndarray
            Response function kernel, shape (n_quad, n_quad)
        """
        assert isinstance(full_eigenvalues, np.ndarray), \
            FULL_EIGENVALUES_TYPE_ERROR_MESSAGE.format(type(full_eigenvalues))
        assert isinstance(full_orbitals, np.ndarray), \
            FULL_ORBITALS_TYPE_ERROR_MESSAGE.format(type(full_orbitals))
        assert isinstance(full_l_terms, np.ndarray), \
            FULL_L_TERMS_TYPE_ERROR_MESSAGE.format(type(full_l_terms))
        
        # Initialize chi_0_kernel
        chi_0_kernel = np.zeros((self.n_quad, self.n_quad))
        
        # Get l channel indices for all orbitals
        l_max = np.max(self.occ_l_values)
        l_channel_orbital_indices = np.zeros((l_max + 1, self.n_interior), dtype=np.int32)

        for l in range(l_max + 1):
            l_channel_orbital_indices[l, :] = np.argwhere(full_l_terms == l)[:, 0]
        
        # Get occupied orbitals
        occ_orbitals = full_orbitals[:, :len(self.occ_l_values)]
        
        # Compute chi_0_kernel for each occupied orbital
        for idx in range(len(self.occ_l_values)):
            # Get l and n index
            l_value = self.occ_l_values[idx]
            n_value = self.occ_n_values[idx] - l_value - 1
            
            # Get all orbitals with indices of the same l value (occupied + unoccupied)
            all_orbitals_in_l_channel = full_orbitals[:, l_channel_orbital_indices[l_value, :]]
            
            # Get the difference of eigenvalues
            l_channel_eigenvalues = full_eigenvalues[l_channel_orbital_indices[l_value, :]]
            diff_eigenvalues = l_channel_eigenvalues.reshape(-1, 1) - l_channel_eigenvalues.reshape(1, -1)
            #   Only compute where |ε_i - ε_j| > threshold, set to 0 otherwise (handles diagonal and degenerate cases)
            threshold = 1e-12
            one_over_diff_eigenvalues = np.divide(1.0, diff_eigenvalues, 
                                                  out=np.zeros_like(diff_eigenvalues),
                                                  where=np.abs(diff_eigenvalues) > threshold)

            # Get the green function block
            _exchange_green_block = np.einsum('ji,ki,i->jk',
                all_orbitals_in_l_channel,
                all_orbitals_in_l_channel,
                one_over_diff_eigenvalues[n_value, :],
                optimize=True
            )
            
            # Get the orbital (on FE grid)
            orbital = occ_orbitals[:, idx]
            
            # Update chi_0_kernel
            # Use interpolated orbital and Green function block (both on quadrature grid)
            chi_0_kernel += 2 * self.occupations[idx] * np.einsum('k,kj,j->kj',
                orbital,
                _exchange_green_block,
                orbital,
                optimize=True
            )
        
        return chi_0_kernel


    def compute_dielectric_matrix(
        self,
        full_eigenvalues : np.ndarray,
        full_orbitals    : np.ndarray,
        full_l_terms     : np.ndarray,
    ) -> np.ndarray:
        """
        Compute dielectric matrix for preconditioning.
        
        The dielectric matrix is defined as:
            ε(r, r') = δ(r - r') - (1/r²) χ₀(r, r') @ ν(r, r')
        
        where ν is the Coulomb kernel.
        
        Parameters
        ----------
        full_eigenvalues : np.ndarray
            Full eigenvalues (occupied + unoccupied), shape (n_states,)
        full_orbitals : np.ndarray
            Full orbitals at quadrature nodes, shape (n_states, n_quad)
        full_l_terms : np.ndarray
            Angular momentum quantum numbers, shape (n_states,)
        
        Returns
        -------
        dielectric_matrix : np.ndarray
            Dielectric matrix, shape (n_quad, n_quad)
        """
        
        # Compute chi_0_kernel
        chi_0_kernel = self.compute_chi_0_kernel(
            full_eigenvalues,
            full_orbitals,
            full_l_terms
        )
        
        # Compute Coulomb kernel ν(r, r') = 1 / max(r, r')
        r_quad = self.quadrature_nodes
        nu = 1.0 / np.maximum(r_quad.reshape(-1, 1), r_quad.reshape(1, -1))
        
        # Compute dielectric matrix
        # Match reference: np.eye(N_q) - (1/r_quad.reshape(-1,1)**2)*(chi_0_kernel*w)@nu*r_quad**2*w
        dielectric_matrix = np.eye(self.n_quad) - (1.0 / r_quad.reshape(-1, 1)**2) * ((chi_0_kernel * self.quadrature_weights) @ nu) * (r_quad**2 * self.quadrature_weights)


        return dielectric_matrix


