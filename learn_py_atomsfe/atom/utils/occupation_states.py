


import numpy as np
from typing import Optional, Tuple


'''
Input: Atomic Number Z
Output: Occ matrix containing "n" quantum number list in first row, "l"
quantum number in the second row, the corresponding spin-up occupation 
in the third row and spin-down occupation in the fourth row.
'''


# Error messages 
Z_NOT_INT_ERROR = \
    "parameter 'Z' must be an integer, get type {} instead."
Z_NOT_IN_VALID_RANGE_ERROR = \
    "parameter 'Z' must be between 1 and 92 (1-92), get {} instead."
Z_NUCLEAR_NOT_INT_OR_FLOAT_ERROR = \
    "parameter 'z_nuclear' must be an integer or float, get type {} instead."
Z_NUCLEAR_NOT_GREATER_THAN_0_OR_LESS_THAN_93_ERROR = \
    "parameter 'z_nuclear' must be greater than 0 and less than 93 (1-92), get {} instead."

Z_VALENCE_NOT_INT_OR_FLOAT_ERROR = \
    "parameter 'z_valence' must be an integer or float, get type {} instead."
Z_VALENCE_NOT_GREATER_THAN_0_OR_LESS_THAN_93_ERROR = \
    "parameter 'z_valence' must be greater than 0 and less than 93 (1-92), get {} instead."
ALL_ELECTRON_FLAG_NOT_BOOL_ERROR = \
    "parameter 'all_electron_flag' must be a boolean, get type {} instead."    
N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR = \
    "parameter 'n_electrons' must be an integer or float, get type {} instead."
N_ELECTRONS_NOT_GREATER_THAN_0_ERROR = \
    "parameter 'n_electrons' must be greater than 0, get {} instead."
N_ELECTRONS_NOT_LESS_THAN_OR_EQUAL_TO_92_ERROR = \
    "parameter 'n_electrons' must be less than or equal to 92, get {} instead."

CHARGE_SYSTEMS_NOT_SUPPORTED_FOR_PSEUDOPOTENTIAL_CALCULATION_ERROR = \
    "Charged systems are not supported with pseudopotentials. Use all-electron calculations for non-neutral systems."
Z_NUCLEAR_NOT_INTEGER_VALUED_FOR_PSEUDOPOTENTIAL_CALCULATION_ERROR = \
    "parameter 'z_nuclear' must be integer-valued for pseudopotential calculations, get {} instead."
TOTAL_OCCUPATION_NUMBERS_DO_NOT_MATCH_THE_NUMBER_OF_ELECTRONS_ERROR = \
    "Total occupation numbers do not match the number of electrons {} != {}, this should not happen."



HARDCORED_OCCUPATION_EXCEPTION_ORBITAL_INDEX_DICT = {
    # The key is the atomic number, the value is the index of the orbital 
    #   that should be changed when the number of electrons is fractional.
    24: -2,  # Cr, 3d
    29: -2,  # Cu, 3d
    41: -2,  # Nb, 4d
    44: -2,  # Ru, 4d
    46: -2,  # Pd, 4p
    58: -5,  # Ce, 4f
    59: -5,  # Pr, 4d
    64: -2,  # Gd, 5d
    65: -5,  # Tb, 4d
    78: -2,  # Pt, 5d
    91: -5,  # Pa, 5f
}



def get_neutral_occupation_states(Z : int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Get occupation states for given atomic number when the system is electrically neutral.

    Parameters
    ----------
    Z : int
        Atomic number

    Returns
    -------
    n_quantum : np.ndarray
        Principal quantum number n for each orbital
    l_quantum : np.ndarray
        Angular momentum quantum number l for each orbital
    s_quantum_up : np.ndarray
        Spin-up occupation for each orbital
    s_quantum_down : np.ndarray
        Spin-down occupation for each orbital
    """
    # Type checking
    assert isinstance(Z, int), \
        Z_NOT_INT_ERROR.format(type(Z))
  
    if Z == 1:
        n_quantum = np.array([1])
        l_quantum = np.array([0])
        s_quantum_up = np.array([1])
        s_quantum_down = np.array([0])
    elif Z == 2:
        n_quantum = np.array([1])
        l_quantum = np.array([0])
        s_quantum_up = np.array([1])
        s_quantum_down = np.array([1])
    elif Z == 3:
        n_quantum = np.array([1,2])
        l_quantum = np.array([0,0])
        s_quantum_up = np.array([1,1])
        s_quantum_down = np.array([1,0])
    elif Z == 4:
        n_quantum = np.array([1,2])
        l_quantum = np.array([0,0])
        s_quantum_up = np.array([1,1])
        s_quantum_down = np.array([1,1])
    elif Z == 5:
        n_quantum = np.array([1,2,2])
        l_quantum = np.array([0,0,1])
        s_quantum_up = np.array([1,1,1])
        s_quantum_down = np.array([1,1,0])
    elif Z == 6:
        n_quantum = np.array([1,2,2])
        l_quantum = np.array([0,0,1])
        s_quantum_up = np.array([1,1,2])
        s_quantum_down = np.array([1,1,0])
    elif Z == 7:
        n_quantum = np.array([1,2,2])
        l_quantum = np.array([0,0,1])
        s_quantum_up = np.array([1,1,3])
        s_quantum_down = np.array([1,1,0])
    elif Z == 8:
        n_quantum = np.array([1,2,2])
        l_quantum = np.array([0,0,1])
        s_quantum_up = np.array([1,1,3])
        s_quantum_down = np.array([1,1,1])
    elif Z == 9:
        n_quantum = np.array([1,2,2])
        l_quantum = np.array([0,0,1])
        s_quantum_up = np.array([1,1,3])
        s_quantum_down = np.array([1,1,2])
    elif Z == 10:
        n_quantum = np.array([1,2,2])
        l_quantum = np.array([0,0,1])
        s_quantum_up = np.array([1,1,3])
        s_quantum_down = np.array([1,1,3])
    elif Z == 11:
        n_quantum = np.array([1,2,2,3])
        l_quantum = np.array([0,0,1,0])
        s_quantum_up = np.array([1,1,3,1])
        s_quantum_down = np.array([1,1,3,0])
    elif Z == 12:
        n_quantum = np.array([1,2,2,3])
        l_quantum = np.array([0,0,1,0])
        s_quantum_up = np.array([1,1,3,1])
        s_quantum_down = np.array([1,1,3,1])
    elif Z == 13:
        n_quantum = np.array([1,2,2,3,3])
        l_quantum = np.array([0,0,1,0,1])
        s_quantum_up = np.array([1,1,3,1,1])
        s_quantum_down = np.array([1,1,3,1,0])
    elif Z == 14:
        n_quantum = np.array([1,2,2,3,3])
        l_quantum = np.array([0,0,1,0,1])
        s_quantum_up = np.array([1,1,3,1,2])
        s_quantum_down = np.array([1,1,3,1,0])
    elif Z == 15:
        n_quantum = np.array([1,2,2,3,3])
        l_quantum = np.array([0,0,1,0,1])
        s_quantum_up = np.array([1,1,3,1,3])
        s_quantum_down = np.array([1,1,3,1,0])
    elif Z == 16:
        n_quantum = np.array([1,2,2,3,3])
        l_quantum = np.array([0,0,1,0,1])
        s_quantum_up = np.array([1,1,3,1,3])
        s_quantum_down = np.array([1,1,3,1,1])
    elif Z == 17:
        n_quantum = np.array([1,2,2,3,3])
        l_quantum = np.array([0,0,1,0,1])
        s_quantum_up = np.array([1,1,3,1,3])
        s_quantum_down = np.array([1,1,3,1,2])
    elif Z == 18:
        n_quantum = np.array([1,2,2,3,3])
        l_quantum = np.array([0,0,1,0,1])
        s_quantum_up = np.array([1,1,3,1,3])
        s_quantum_down = np.array([1,1,3,1,3])
    elif Z == 19:
        n_quantum = np.array([1,2,2,3,3,4])
        l_quantum = np.array([0,0,1,0,1,0])
        s_quantum_up = np.array([1,1,3,1,3,1])
        s_quantum_down = np.array([1,1,3,1,3,0])
    elif Z == 20:
        n_quantum = np.array([1,2,2,3,3,4])
        l_quantum = np.array([0,0,1,0,1,0])
        s_quantum_up = np.array([1,1,3,1,3,1])
        s_quantum_down = np.array([1,1,3,1,3,1])
    elif Z == 21:
        n_quantum = np.array([1, 2, 2, 3, 3, 3, 4])
        l_quantum = np.array([0, 0, 1, 0, 1, 2, 0])
        s_quantum_up = np.array([1, 1, 3, 1, 3, 1, 1])
        s_quantum_down = np.array([1, 1, 3, 1, 3, 0, 1])
    elif Z == 22:
        n_quantum = np.array([1, 2, 2, 3, 3, 3, 4])
        l_quantum = np.array([0, 0, 1, 0, 1, 2, 0])
        s_quantum_up = np.array([1, 1, 3, 1, 3, 2, 1])
        s_quantum_down = np.array([1, 1, 3, 1, 3, 0, 1])
    elif Z == 23:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 0, 1 ])
    elif Z == 24:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 0, 0 ])
    elif Z == 25:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 0, 1 ])
    elif Z == 26:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 1, 1 ])
    elif Z == 27:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 2, 1 ])
    elif Z == 28:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 3, 1 ])
    elif Z == 29:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 0 ])
    elif Z == 30:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1 ])
    elif Z == 31:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 0 ])
    elif Z == 32:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 2 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 0 ])
    elif Z == 33:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 0 ])
    elif Z == 34:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 1 ])
    elif Z == 35:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 2 ])
    elif Z == 36:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3 ])
    elif Z == 37:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 0 ])
    elif Z == 38:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 1 ])
    elif Z == 39:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 0, 1 ])
    elif Z == 40:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 2, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 0, 1 ])
    elif Z == 41:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 4, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 0, 0 ])
    elif Z == 42:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 0, 0 ])
    elif Z == 43:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 0, 1 ])
    elif Z == 44:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 2, 0 ])
    elif Z == 45:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 3, 0 ])
    elif Z == 46:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5 ])
    elif Z == 47:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0 ])
    elif Z == 48:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1 ])
    elif Z == 49:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 0 ])
    elif Z == 50:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 2 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 0 ])
    elif Z == 51:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 0 ])
    elif Z == 52:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 1 ])
    elif Z == 53:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 2 ])
    elif Z == 54:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3 ])
    elif Z == 55:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3, 0 ])
    elif Z == 56:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3, 1 ])
    elif Z == 57:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 3, 0, 1 ])
    elif Z == 58:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 1, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 0, 1 ])
    elif Z == 59:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 3, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 1 ])
    elif Z == 60:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 4, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 1 ])
    elif Z == 61:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 1 ])
    elif Z == 62:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 6, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 1 ])
    elif Z == 63:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 1 ])
    elif Z == 64:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 0, 1, 3, 0, 1 ])
    elif Z == 65:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 2, 1, 3, 1 ])
    elif Z == 66:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 3, 1, 3, 1 ])
    elif Z == 67:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 4, 1, 3, 1 ])
    elif Z == 68:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 5, 1, 3, 1 ])
    elif Z == 69:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 6, 1, 3, 1 ])
    elif Z == 70:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1 ])
    elif Z == 71:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 0, 1 ])
    elif Z == 72:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 2, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 0, 1 ])
    elif Z == 73:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 0, 1 ])
    elif Z == 74:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 4, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 0, 1 ])
    elif Z == 75:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 0, 1 ])
    elif Z == 76:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 1, 1 ])
    elif Z == 77:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 2, 1 ])
    elif Z == 78:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 4, 0 ])
    elif Z == 79:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 0 ])
    elif Z == 80:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1 ])
    elif Z == 81:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 0 ])
    elif Z == 82:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 2 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 0 ])
    elif Z == 83:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 0 ])
    elif Z == 84:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 1 ])
    elif Z == 85:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 2 ])
    elif Z == 86:
        n_quantum = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6 ])
        l_quantum = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1 ])
        s_quantum_up = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3 ])
    elif Z == 87:
        n_quantum      = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 ])
        l_quantum      = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 ])
        s_quantum_up   = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 0 ])
    elif Z == 88:
        n_quantum      = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 ])
        l_quantum      = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 ])
        s_quantum_up   = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 1 ])
    elif Z == 89:
        n_quantum      = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7 ])
        l_quantum      = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up   = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 0, 1 ])
    elif Z == 90:
        n_quantum      = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7 ])
        l_quantum      = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0 ])
        s_quantum_up   = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 2, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 1, 3, 0, 1 ])
    elif Z == 91:
        n_quantum      = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7 ])
        l_quantum      = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up   = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 2, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 0, 1, 3, 0, 1 ])
    elif Z == 92:
        n_quantum      = np.array([ 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7 ])
        l_quantum      = np.array([ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0 ])
        s_quantum_up   = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 3, 1, 3, 1, 1 ])
        s_quantum_down = np.array([ 1, 1, 3, 1, 3, 5, 1, 3, 5, 7, 1, 3, 5, 0, 1, 3, 0, 1 ])
    else:
        raise ValueError(Z_NOT_IN_VALID_RANGE_ERROR.format(Z))
    
    return n_quantum, l_quantum, s_quantum_up, s_quantum_down 



def get_fraction_occupation_states(n_electrons : float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Get occupation states for fractional number of electrons.

    Parameters
    ----------
    n_electrons : float
        Number of electrons in the system

    Returns
    -------
    n_quantum : np.ndarray[int]
        Principal quantum number n for each orbital
    l_quantum : np.ndarray[int]
        Angular momentum quantum number l for each orbital
    s_quantum_up : np.ndarray[float]
        Spin-up occupation for each orbital
    s_quantum_down : np.ndarray[float]
        Spin-down occupation for each orbital
    """

    # Type checking
    assert isinstance(n_electrons, int | float), \
        N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(type(n_electrons))
    assert n_electrons > 0, \
        N_ELECTRONS_NOT_GREATER_THAN_0_ERROR.format(n_electrons)
    assert n_electrons <= 92, \
        N_ELECTRONS_NOT_LESS_THAN_OR_EQUAL_TO_92_ERROR.format(n_electrons)

    # Use the smallest integer >= n_electrons as the reference configuration.
    ceil_electrons = int(np.ceil(n_electrons))
    fractional_delta = n_electrons - ceil_electrons  # <= 0 for fractional (removing charge)

    n_quantum, l_quantum, occ_spin_up, occ_spin_down = get_neutral_occupation_states(ceil_electrons)
    occ_spin_up   = occ_spin_up.astype(float)
    occ_spin_down = occ_spin_down.astype(float)

    if n_electrons.is_integer():
        assert fractional_delta == 0
    elif ceil_electrons not in HARDCORED_OCCUPATION_EXCEPTION_ORBITAL_INDEX_DICT:
        if ceil_electrons == 1:
            # Hydrogen-like case: only one orbital exists.
            occ_spin_up   = np.array([n_electrons], dtype=float)
            occ_spin_down = np.array([0.0], dtype=float)
        else:
            # Align Z-1 orbitals to the Z orbital list, then locate the changed orbital.
            (n_quantum_prev, l_quantum_prev, occ_spin_up_prev, occ_spin_down_prev) = get_neutral_occupation_states(ceil_electrons - 1)

            occ_spin_up_prev_aligned   = np.zeros_like(occ_spin_up, dtype=float)
            occ_spin_down_prev_aligned = np.zeros_like(occ_spin_down, dtype=float)

            prev_idx = 0
            curr_idx = 0
            # Two-pointer walk to align Z-1 orbitals into the Z orbital list.
            # When (n,l) matches, copy the previous occupations into the current index.
            while prev_idx < len(n_quantum_prev) and curr_idx < len(n_quantum):
                if (n_quantum_prev[prev_idx] == n_quantum[curr_idx]) and \
                   (l_quantum_prev[prev_idx] == l_quantum[curr_idx]):
                    occ_spin_up_prev_aligned[curr_idx] = occ_spin_up_prev[prev_idx]
                    occ_spin_down_prev_aligned[curr_idx] = occ_spin_down_prev[prev_idx]
                    prev_idx += 1
                    curr_idx += 1
                else:
                    # The current Z list has an extra (new) orbital not in Z-1.
                    # Leave its aligned occupation as zero and advance only curr_idx.
                    curr_idx += 1

            # Calculate the difference between the current and previous occupations.
            diff_up   = occ_spin_up - occ_spin_up_prev_aligned
            diff_down = occ_spin_down - occ_spin_down_prev_aligned

            # Find the indices of the orbitals that have changed.
            changed_up   = np.argwhere(diff_up != 0)[:, 0]
            changed_down = np.argwhere(diff_down != 0)[:, 0]

            if len(changed_up) > 0:
                occ_spin_up[changed_up] = occ_spin_up[changed_up] + fractional_delta
            elif len(changed_down) > 0:
                occ_spin_down[changed_down] = occ_spin_down[changed_down] + fractional_delta

            if len(changed_up) > 0 and len(changed_down) > 0:
                raise ValueError("Ambiguous fractional occupation change in both spins, this should not happen.")
    else:
        # Hardcoded exceptions for discontinuous occupation changes.
        assert ceil_electrons in HARDCORED_OCCUPATION_EXCEPTION_ORBITAL_INDEX_DICT, \
            "Invalid atomic number {} for fractional occupation states, this should not happen.".format(ceil_electrons)
        offset = HARDCORED_OCCUPATION_EXCEPTION_ORBITAL_INDEX_DICT[ceil_electrons]
        occ_spin_up[offset] = occ_spin_up[offset] + fractional_delta

    return n_quantum, l_quantum, occ_spin_up, occ_spin_down


class OccupationInfo:
    """
    Occupation information for atomic states.
    """
    z_valence                  : int | float  # Valence charge (for pseudopotential)
    z_nuclear                  : int | float  # True nuclear charge of the atom
    n_electrons                : int | float  # Number of electrons in the system, can be fractional
    all_electron_flag          : bool         # Whether to use all-electron or pseudopotential
    occ_n                      : np.ndarray   # Principal quantum number n for each orbital
    occ_l                      : np.ndarray   # Angular momentum quantum number l for each orbital
    occ_spin_up                : np.ndarray   # Spin-up occupation for each orbital
    occ_spin_down              : np.ndarray   # Spin-down occupation for each orbital
    occ_spin_up_plus_spin_down : np.ndarray   # Total occupation (spin-up + spin-down)


    def __init__(self, 
        z_nuclear         : int | float,                  # True nuclear charge (atomic number)
        z_valence         : int | float,                  # Valence charge (for pseudopotential Coulomb tail)
        all_electron_flag : bool,                         # Whether to use all-electron or pseudopotential
        n_electrons       : Optional[int | float] = None, # Number of electrons in the system, can be fractional
    ):
        """
        Initialize occupation information.
        
        Parameters
        ----------
        z_nuclear : int | float
            True nuclear charge of the atom (atomic number)
        z_valence : int | float
            Valence charge for pseudopotential calculations
        all_electron_flag : bool
            True for all-electron, False for pseudopotential
        n_electrons : int | float, optional
            Number of electrons in the system, can be fractional, by default, set to float(z_nuclear)
        """

        self.z_nuclear         : float = z_nuclear
        self.z_valence         : float = z_valence
        self.n_electrons       : float = n_electrons
        self.all_electron_flag : bool  = all_electron_flag
        self._set_and_check_initial_parameters()


        if all_electron_flag:

            n_quantum, l_quantum, s_quantum_up, s_quantum_down = get_fraction_occupation_states(self.n_electrons)
            self.occ_n = n_quantum
            self.occ_l = l_quantum
            self.occ_spin_up   = s_quantum_up
            self.occ_spin_down = s_quantum_down
        else:
            # For pseudopotential: check if z_nuclear is integer-valued and equal to n_electrons
            assert self.z_nuclear.is_integer(), \
                Z_NUCLEAR_NOT_INTEGER_VALUED_FOR_PSEUDOPOTENTIAL_CALCULATION_ERROR.format(self.z_nuclear)
            assert self.z_nuclear == self.n_electrons, \
                CHARGE_SYSTEMS_NOT_SUPPORTED_FOR_PSEUDOPOTENTIAL_CALCULATION_ERROR.format(self.n_electrons)

            n_quantum, l_quantum, s_quantum_up, s_quantum_down = get_neutral_occupation_states(int(self.z_nuclear))

            # select only valence electrons
            n_core_electrons = z_nuclear - z_valence
            orbital_occupation_numbers = s_quantum_up + s_quantum_down
            cumulative_occupation = np.cumsum(orbital_occupation_numbers)
            valence_orbitals_indices = np.where(cumulative_occupation > n_core_electrons)[0]
            self.occ_n = n_quantum[valence_orbitals_indices]
            self.occ_l = l_quantum[valence_orbitals_indices]
            self.occ_spin_up   = s_quantum_up[valence_orbitals_indices]
            self.occ_spin_down = s_quantum_down[valence_orbitals_indices]
            
        self.occ_spin_up_plus_spin_down = self.occ_spin_up + self.occ_spin_down

        # Check if the total occupation numbers match the number of electrons
        if self.all_electron_flag:
            assert np.sum(self.occ_spin_up_plus_spin_down) == self.n_electrons, \
                TOTAL_OCCUPATION_NUMBERS_DO_NOT_MATCH_THE_NUMBER_OF_ELECTRONS_ERROR.format(np.sum(self.occ_spin_up_plus_spin_down), self.n_electrons)
        else:
            assert np.sum(self.occ_spin_up_plus_spin_down) == self.z_valence, \
                TOTAL_OCCUPATION_NUMBERS_DO_NOT_MATCH_THE_NUMBER_OF_ELECTRONS_ERROR.format(np.sum(self.occ_spin_up_plus_spin_down), self.z_valence)



    def _set_and_check_initial_parameters(self):
        """
        Set and check initial parameters.
        """
        # Type checking and conversion
        # z_nuclear
        assert isinstance(self.z_nuclear, (int, float)), \
            Z_NUCLEAR_NOT_INT_OR_FLOAT_ERROR.format(type(self.z_nuclear))
        assert self.z_nuclear > 0 and self.z_nuclear < 93, \
            Z_NUCLEAR_NOT_GREATER_THAN_0_OR_LESS_THAN_93_ERROR.format(self.z_nuclear)
        self.z_nuclear = float(self.z_nuclear)

        # z_valence
        assert isinstance(self.z_valence, (int, float)), \
            Z_VALENCE_NOT_INT_OR_FLOAT_ERROR.format(type(self.z_valence))
        assert self.z_valence > 0 and self.z_valence < 93, \
            Z_VALENCE_NOT_GREATER_THAN_0_OR_LESS_THAN_93_ERROR.format(self.z_valence)
        self.z_valence = float(self.z_valence)

        # all electron flag
        assert isinstance(self.all_electron_flag, bool), \
            ALL_ELECTRON_FLAG_NOT_BOOL_ERROR.format(type(self.all_electron_flag))
        
        # n_electrons
        if self.n_electrons is None:
            self.n_electrons = self.z_nuclear
        else:
            assert isinstance(self.n_electrons, (int, float)), \
                N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(type(self.n_electrons))
            assert self.n_electrons > 0, \
                N_ELECTRONS_NOT_GREATER_THAN_0_ERROR.format(self.n_electrons)
            self.n_electrons = float(self.n_electrons)


    @property
    def n_free_electrons(self) -> float:
        """
        Total number of free electrons in the system.

        This is the sum of spin-up and spin-down occupations over the active
        orbital set:
        - All-electron calculations: sum over all-electron occupations.
        - Pseudopotential calculations: sum over valence occupations.
        """
        return np.sum(self.occ_spin_up_plus_spin_down)

    @property
    def n_free_electrons_up(self) -> float:
        """
        Number of spin-up free electrons in the system.

        This is the sum of spin-up occupations over the active orbital set:
        - All-electron calculations: sum over all-electron spin-up occupations.
        - Pseudopotential calculations: sum over valence spin-up occupations.
        """
        return np.sum(self.occ_spin_up)
    
    @property
    def n_free_electrons_dn(self) -> float:
        """
        Number of spin-down free electrons in the system.

        This is the sum of spin-down occupations over the active orbital set:
        - All-electron calculations: sum over all-electron spin-down occupations.
        - Pseudopotential calculations: sum over valence spin-down occupations.
        """
        return np.sum(self.occ_spin_down)


    @property
    def occupations(self) -> np.ndarray:
        """
        Total occupation numbers (spin-up + spin-down) for each orbital.
        Alias for occ_spin_up_plus_spin_down for cleaner API.
        """
        return self.occ_spin_up_plus_spin_down
    
    @property
    def occupations_up(self) -> np.ndarray:
        """
        Spin-up occupation numbers for each orbital.
        """
        return self.occ_spin_up
    
    @property
    def occupations_dn(self) -> np.ndarray:
        """
        Spin-down occupation numbers for each orbital.
        """
        return self.occ_spin_down




    @property
    def l_values(self) -> np.ndarray:
        """
        Angular momentum quantum numbers for each orbital.
        Alias for occ_l for cleaner API.
        """
        return self.occ_l
    
    @property
    def n_values(self) -> np.ndarray:
        """
        Principal quantum numbers for each orbital.
        Alias for occ_n for cleaner API.
        """
        return self.occ_n
    
    @property
    def unique_l_values(self) -> np.ndarray:
        """Get unique angular momentum quantum numbers present in occupied states."""
        return np.unique(self.occ_l)
    

    @property
    def n_states(self) -> int:
        """Get total number of occupied states."""
        return len(self.occ_n)


    def n_states_for_l(self, l: int) -> int:
        """
        Get number of occupied states for a given angular momentum quantum number.
        
        Parameters
        ----------
        l : int
            Angular momentum quantum number
        
        Returns
        -------
        n_states : int
            Number of states with this l value
        """
        return np.sum(self.occ_l == l)



    @property
    def closed_shell_flag(self) -> bool:
        """
        Check if the atom is closed-shell.
        """
        for n, l, spin_up, spin_down in zip(self.occ_n, self.occ_l, self.occ_spin_up, self.occ_spin_down):
            if spin_up != l * 2 + 1 or spin_down != l * 2 + 1:
                return False
        return True



    def print_info(self):
        def _format_array_for_print(values: np.ndarray) -> np.ndarray:
            """
            Format array for display only.
            If all entries are integer-valued, print as ints.
            """
            array_values = np.asarray(values)
            if np.issubdtype(array_values.dtype, np.integer):
                return array_values
            if np.allclose(array_values, np.rint(array_values), rtol=0.0, atol=1e-12):
                return np.rint(array_values).astype(int)
            return array_values

        print("===========================================================================")
        print("                        OCCUPATION INFORMATION                             ")
        print("===========================================================================")
        print(f"\t z_valence (valence charge) : {self.z_valence}")
        print(f"\t z_nuclear (nuclear charge) : {self.z_nuclear}")
        print(f"\t n_electrons                : {self.n_electrons}")
        print(f"\t all_electron_flag          : {self.all_electron_flag}")
        print(f"\t occ_n                      : {_format_array_for_print(self.occ_n)}")
        print(f"\t occ_l                      : {_format_array_for_print(self.occ_l)}")
        print(f"\t occ_spin_up                : {_format_array_for_print(self.occ_spin_up)}")
        print(f"\t occ_spin_down              : {_format_array_for_print(self.occ_spin_down)}")
        print(f"\t occ_spin_up_plus_spin_down : {_format_array_for_print(self.occ_spin_up_plus_spin_down)}")
        print()




if __name__ == "__main__":
    for atomic_number in range(1, 93):
        occupation_info = OccupationInfo(z_nuclear=atomic_number, z_valence=atomic_number, all_electron_flag=True)
        if occupation_info.closed_shell_flag:
            print(f"atomic_number = {atomic_number} is closed shell")
        else:
            pass



