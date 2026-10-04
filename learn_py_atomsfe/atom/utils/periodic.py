
# error messages
ATOMIC_NUMBER_NOT_INTEGER_ERROR = \
    "Atomic number must be an integer, get type {} instead."
ATOMIC_NUMBER_NOT_SUPPORTED_ERROR = \
    "Atomic number {} is not supported. Please use an integer between 1 and 93."
ATOM_NAME_NOT_STRING_ERROR = \
    "Atom name must be a string, get type {} instead."
ATOM_NAME_NOT_SUPPORTED_ERROR = \
    "Atom named '{}' is not supported."


# Common oxidation states (valences) for each element
# Values are lists of common oxidation states, ordered by frequency/importance
# Note: For charged systems, n_electrons = z_nuclear - oxidation_state
COMMON_OXIDATION_STATES = {
    1:  [1, -1, 0],               # H: +1 (most common), -1 (hydride), 0 (neutral)
    2:  [0],                      # He: 0 (noble gas, rarely forms compounds)
    3:  [1, 0],                   # Li: +1, 0
    4:  [2, 0],                   # Be: +2, 0
    5:  [3, 0],                   # B: +3, 0
    6:  [4, 2, -4, 0],            # C: +4, +2, -4, 0
    7:  [3, 5, -3, 0],            # N: +3, +5, -3, 0
    8:  [-2, 0],                  # O: -2, 0
    9:  [-1, 0],                  # F: -1, 0
    10: [0],                      # Ne: 0
    11: [1, 0],                   # Na: +1, 0
    12: [2, 0],                   # Mg: +2, 0
    13: [3, 0],                   # Al: +3, 0
    14: [4, 2, -4, 0],            # Si: +4, +2, -4, 0
    15: [5, 3, -3, 0],            # P: +5, +3, -3, 0
    16: [6, 4, 2, -2, 0],         # S: +6, +4, +2, -2, 0
    17: [7, 5, 3, 1, -1, 0],      # Cl: +7, +5, +3, +1, -1, 0
    18: [0],                      # Ar: 0
    19: [1, 0],                   # K: +1, 0
    20: [2, 0],                   # Ca: +2, 0
    21: [3, 2, 0],                # Sc: +3, +2, 0
    22: [4, 3, 2, 0],             # Ti: +4, +3, +2, 0
    23: [5, 4, 3, 2, 0],          # V: +5, +4, +3, +2, 0
    24: [6, 3, 2, 0],             # Cr: +6, +3, +2, 0
    25: [7, 6, 4, 3, 2, 0],       # Mn: +7, +6, +4, +3, +2, 0
    26: [3, 2, 0],                # Fe: +3, +2, 0
    27: [3, 2, 0],                # Co: +3, +2, 0
    28: [2, 3, 0],                # Ni: +2, +3, 0
    29: [2, 1, 0],                # Cu: +2, +1, 0
    30: [2, 0],                   # Zn: +2, 0
    31: [3, 1, 0],                # Ga: +3, +1, 0
    32: [4, 2, 0],                # Ge: +4, +2, 0
    33: [5, 3, -3, 0],            # As: +5, +3, -3, 0
    34: [6, 4, 2, -2, 0],         # Se: +6, +4, +2, -2, 0
    35: [7, 5, 3, 1, -1, 0],      # Br: +7, +5, +3, +1, -1, 0
    36: [0, 2],                   # Kr: 0, +2 (rare)
    37: [1, 0],                   # Rb: +1, 0
    38: [2, 0],                   # Sr: +2, 0
    39: [3, 0],                   # Y: +3, 0
    40: [4, 0],                   # Zr: +4, 0
    41: [5, 4, 3, 0],             # Nb: +5, +4, +3, 0
    42: [6, 5, 4, 3, 0],          # Mo: +6, +5, +4, +3, 0
    43: [7, 6, 4, 0],             # Tc: +7, +6, +4, 0
    44: [4, 3, 2, 0],             # Ru: +4, +3, +2, 0
    45: [3, 4, 0],                # Rh: +3, +4, 0
    46: [2, 4, 0],                # Pd: +2, +4, 0
    47: [1, 2, 0],                # Ag: +1, +2, 0
    48: [2, 0],                   # Cd: +2, 0
    49: [3, 1, 0],                # In: +3, +1, 0
    50: [4, 2, 0],                # Sn: +4, +2, 0
    51: [5, 3, -3, 0],            # Sb: +5, +3, -3, 0
    52: [6, 4, 2, -2, 0],         # Te: +6, +4, +2, -2, 0
    53: [7, 5, 3, 1, -1, 0],      # I: +7, +5, +3, +1, -1, 0
    54: [0, 2, 4, 6, 8],          # Xe: 0, +2, +4, +6, +8 (rare)
    55: [1, 0],                   # Cs: +1, 0
    56: [2, 0],                   # Ba: +2, 0
    57: [3, 0],                   # La: +3, 0
    58: [4, 3, 0],                # Ce: +4, +3, 0
    59: [4, 3, 0],                # Pr: +4, +3, 0
    60: [3, 0],                   # Nd: +3, 0
    61: [3, 0],                   # Pm: +3, 0
    62: [3, 2, 0],                # Sm: +3, +2, 0
    63: [3, 2, 0],                # Eu: +3, +2, 0
    64: [3, 0],                   # Gd: +3, 0
    65: [4, 3, 0],                # Tb: +4, +3, 0
    66: [3, 0],                   # Dy: +3, 0
    67: [3, 0],                   # Ho: +3, 0
    68: [3, 0],                   # Er: +3, 0
    69: [3, 2, 0],                # Tm: +3, +2, 0
    70: [3, 2, 0],                # Yb: +3, +2, 0
    71: [3, 0],                   # Lu: +3, 0
    72: [4, 0],                   # Hf: +4, 0
    73: [5, 4, 0],                # Ta: +5, +4, 0
    74: [6, 5, 4, 0],             # W: +6, +5, +4, 0
    75: [7, 6, 5, 4, 0],          # Re: +7, +6, +5, +4, 0
    76: [8, 6, 4, 3, 0],          # Os: +8, +6, +4, +3, 0
    77: [4, 3, 0],                # Ir: +4, +3, 0
    78: [4, 2, 0],                # Pt: +4, +2, 0
    79: [3, 1, 0],                # Au: +3, +1, 0
    80: [2, 1, 0],                # Hg: +2, +1, 0
    81: [3, 1, 0],                # Tl: +3, +1, 0
    82: [4, 2, 0],                # Pb: +4, +2, 0
    83: [5, 3, 0],                # Bi: +5, +3, 0
    84: [4, 2, -2, 0],            # Po: +4, +2, -2, 0
    85: [5, 3, 1, -1, 0],         # At: +5, +3, +1, -1, 0
    86: [0, 2, 4, 6],             # Rn: 0, +2, +4, +6 (rare)
    87: [1, 0],                   # Fr: +1, 0
    88: [2, 0],                   # Ra: +2, 0
    89: [3, 0],                   # Ac: +3, 0
    90: [4, 3, 0],                # Th: +4, +3, 0
    91: [5, 4, 3, 0],             # Pa: +5, +4, +3, 0
    92: [6, 5, 4, 3, 0],          # U: +6, +5, +4, +3, 0
}


def get_default_charged_dataset_atomic_number_and_n_electrons_list():
    default_charged_dataset_atomic_number_list = []
    default_charged_dataset_n_electrons_list = []
    for atomic_number in range(1, 93):
        for oxidation_state in COMMON_OXIDATION_STATES[atomic_number]:
            n_electrons = atomic_number - oxidation_state
            if n_electrons > 0:
                default_charged_dataset_atomic_number_list.append(atomic_number)
                default_charged_dataset_n_electrons_list.append(n_electrons)

    return default_charged_dataset_atomic_number_list, default_charged_dataset_n_electrons_list





def atomic_number_to_name(atomic_number: int) -> str:
    assert isinstance(atomic_number, int), \
        ATOMIC_NUMBER_NOT_INTEGER_ERROR.format(type(atomic_number))
    
    if atomic_number == 1: return "H"
    elif atomic_number == 2: return "He"
    elif atomic_number == 3: return "Li"
    elif atomic_number == 4: return "Be"
    elif atomic_number == 5: return "B"
    elif atomic_number == 6: return "C"
    elif atomic_number == 7: return "N"
    elif atomic_number == 8: return "O"
    elif atomic_number == 9: return "F"
    elif atomic_number == 10: return "Ne"
    elif atomic_number == 11: return "Na"
    elif atomic_number == 12: return "Mg"
    elif atomic_number == 13: return "Al"
    elif atomic_number == 14: return "Si"
    elif atomic_number == 15: return "P"
    elif atomic_number == 16: return "S"
    elif atomic_number == 17: return "Cl"
    elif atomic_number == 18: return "Ar"
    elif atomic_number == 19: return "K"
    elif atomic_number == 20: return "Ca"
    elif atomic_number == 21: return "Sc"
    elif atomic_number == 22: return "Ti"
    elif atomic_number == 23: return "V"
    elif atomic_number == 24: return "Cr"
    elif atomic_number == 25: return "Mn"
    elif atomic_number == 26: return "Fe"
    elif atomic_number == 27: return "Co"
    elif atomic_number == 28: return "Ni"
    elif atomic_number == 29: return "Cu"
    elif atomic_number == 30: return "Zn"
    elif atomic_number == 31: return "Ga"
    elif atomic_number == 32: return "Ge"
    elif atomic_number == 33: return "As"
    elif atomic_number == 34: return "Se"
    elif atomic_number == 35: return "Br"
    elif atomic_number == 36: return "Kr"
    elif atomic_number == 37: return "Rb"
    elif atomic_number == 38: return "Sr"
    elif atomic_number == 39: return "Y"
    elif atomic_number == 40: return "Zr"
    elif atomic_number == 41: return "Nb"
    elif atomic_number == 42: return "Mo"
    elif atomic_number == 43: return "Tc"
    elif atomic_number == 44: return "Ru"
    elif atomic_number == 45: return "Rh"
    elif atomic_number == 46: return "Pd"
    elif atomic_number == 47: return "Ag"
    elif atomic_number == 48: return "Cd"
    elif atomic_number == 49: return "In"
    elif atomic_number == 50: return "Sn"
    elif atomic_number == 51: return "Sb"
    elif atomic_number == 52: return "Te"
    elif atomic_number == 53: return "I"
    elif atomic_number == 54: return "Xe"
    elif atomic_number == 55: return "Cs"
    elif atomic_number == 56: return "Ba"
    elif atomic_number == 57: return "La"
    elif atomic_number == 58: return "Ce"
    elif atomic_number == 59: return "Pr"
    elif atomic_number == 60: return "Nd"
    elif atomic_number == 61: return "Pm"
    elif atomic_number == 62: return "Sm"
    elif atomic_number == 63: return "Eu"
    elif atomic_number == 64: return "Gd"
    elif atomic_number == 65: return "Tb"
    elif atomic_number == 66: return "Dy"
    elif atomic_number == 67: return "Ho"
    elif atomic_number == 68: return "Er"
    elif atomic_number == 69: return "Tm"
    elif atomic_number == 70: return "Yb"
    elif atomic_number == 71: return "Lu"
    elif atomic_number == 72: return "Hf"
    elif atomic_number == 73: return "Ta"
    elif atomic_number == 74: return "W"
    elif atomic_number == 75: return "Re"
    elif atomic_number == 76: return "Os"
    elif atomic_number == 77: return "Ir"
    elif atomic_number == 78: return "Pt"
    elif atomic_number == 79: return "Au"
    elif atomic_number == 80: return "Hg"
    elif atomic_number == 81: return "Tl"
    elif atomic_number == 82: return "Pb"
    elif atomic_number == 83: return "Bi"
    elif atomic_number == 84: return "Po"
    elif atomic_number == 85: return "At"
    elif atomic_number == 86: return "Rn"
    elif atomic_number == 87: return "Fr"
    elif atomic_number == 88: return "Ra"
    elif atomic_number == 89: return "Ac"
    elif atomic_number == 90: return "Th"
    elif atomic_number == 91: return "Pa"
    elif atomic_number == 92: return "U"
    elif atomic_number == 93: return "Np"
    else:
        raise ValueError(ATOMIC_NUMBER_NOT_SUPPORTED_ERROR.format(atomic_number))


def name_to_atomic_number(name: str) -> str:
    assert isinstance(name, str), \
        ATOM_NAME_NOT_STRING_ERROR.format(type(name))
    
    if name == "H": return "01"
    elif name == "He": return "02"
    elif name == "Li": return "03"
    elif name == "Be": return "04"
    elif name == "B": return "05"
    elif name == "C": return "06"
    elif name == "N": return "07"
    elif name == "O": return "08"
    elif name == "F": return "09"
    elif name == "Ne": return "10"
    elif name == "Na": return "11"
    elif name == "Mg": return "12"
    elif name == "Al": return "13"
    elif name == "Si": return "14"
    elif name == "P": return "15"
    elif name == "S": return "16"
    elif name == "Cl": return "17"
    elif name == "Ar": return "18"
    elif name == "K": return "19"
    elif name == "Ca": return "20"      
    elif name == "Sc": return "21"
    elif name == "Ti": return "22"
    elif name == "V": return "23"
    elif name == "Cr": return "24"
    elif name == "Mn": return "25"
    elif name == "Fe": return "26"
    elif name == "Co": return "27"
    elif name == "Ni": return "28"
    elif name == "Cu": return "29"
    elif name == "Zn": return "30"
    elif name == "Ga": return "31"
    elif name == "Ge": return "32"
    elif name == "As": return "33"    
    elif name == "Se": return "34"
    elif name == "Br": return "35"
    elif name == "Kr": return "36"
    elif name == "Rb": return "37"
    elif name == "Sr": return "38"
    elif name == "Y": return "39"
    elif name == "Zr": return "40"
    elif name == "Nb": return "41"
    elif name == "Mo": return "42"
    elif name == "Tc": return "43"
    elif name == "Ru": return "44"
    elif name == "Rh": return "45"
    elif name == "Pd": return "46"
    elif name == "Ag": return "47"
    elif name == "Cd": return "48"
    elif name == "In": return "49"
    elif name == "Sn": return "50"
    elif name == "Sb": return "51"
    elif name == "Te": return "52"
    elif name == "I": return "53"
    elif name == "Xe": return "54"
    elif name == "Cs": return "55"
    elif name == "Ba": return "56"
    elif name == "La": return "57"
    elif name == "Ce": return "58"
    elif name == "Pr": return "59"
    elif name == "Nd": return "60"
    elif name == "Pm": return "61"
    elif name == "Sm": return "62"
    elif name == "Eu": return "63"
    elif name == "Gd": return "64"
    elif name == "Tb": return "65"  
    elif name == "Dy": return "66"
    elif name == "Ho": return "67"
    elif name == "Er": return "68"
    elif name == "Tm": return "69"
    elif name == "Yb": return "70"
    elif name == "Lu": return "71"
    elif name == "Hf": return "72"
    elif name == "Ta": return "73"
    elif name == "W": return "74"
    elif name == "Re": return "75"
    elif name == "Os": return "76"
    elif name == "Ir": return "77"  
    elif name == "Pt": return "78"
    elif name == "Au": return "79"
    elif name == "Hg": return "80"
    elif name == "Tl": return "81"
    elif name == "Pb": return "82"
    elif name == "Bi": return "83"
    elif name == "Po": return "84"
    elif name == "At": return "85"
    elif name == "Rn": return "86"
    elif name == "Fr": return "87"
    elif name == "Ra": return "88"
    elif name == "Ac": return "89"
    elif name == "Th": return "90"
    elif name == "Pa": return "91"
    elif name == "U": return "92"
    elif name == "Np": return "93"
    else:
        raise ValueError(ATOMIC_NUMBER_NOT_SUPPORTED_ERROR.format(name))



if __name__ == "__main__":
    print(get_default_charged_dataset_atomic_number_and_n_electrons_list())
    # count the number charged neutal, charged positive, charged negative
    charged_neutal = 0
    charged_positive = 0
    charged_negative = 0
    for atomic_number, n_electrons in zip(get_default_charged_dataset_atomic_number_and_n_electrons_list()[0], get_default_charged_dataset_atomic_number_and_n_electrons_list()[1]):
        if n_electrons- atomic_number == 0:
            charged_neutal += 1
        elif n_electrons- atomic_number > 0:
            charged_positive += 1
        elif n_electrons- atomic_number < 0:
            charged_negative += 1
    print(f"charged_neutal: {charged_neutal}, charged_positive: {charged_positive}, charged_negative: {charged_negative}")