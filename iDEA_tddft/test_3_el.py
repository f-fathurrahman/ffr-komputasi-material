import numpy as np
import matplotlib.pyplot as plt

import iDEA

# A predefined system
#atom = iDEA.system.systems.atom # too large, 4 electrons

#__x2 = np.linspace(-20, 200, 200)
#atom = iDEA.system.System(
#    __x2, -2.0 / (abs(__x2) + 1.0), iDEA.interactions.softened_interaction(__x2), "ud"
#)

# Harmonic potential
#x = np.linspace(-10, 10, 150)
#v_ext = 0.5 * 0.25**2 * x**2
#v_int = iDEA.interactions.softened_interaction(x)
#atom = iDEA.system.System(x, v_ext, v_int, electrons='ud') # also test ud

ω = 1.0
x = np.linspace(-8, 8, 20)
v_ext = 0.5 * ω**2 * x**2
v_int = iDEA.interactions.softened_interaction_alternative(x)
atom = iDEA.system.System(x, v_ext, v_int, electrons='udu', stencil=3)

ground_state = iDEA.methods.interacting.solve(atom, k=0)
print("Total energy = ", ground_state.energy)

