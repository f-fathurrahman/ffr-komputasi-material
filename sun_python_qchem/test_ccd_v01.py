from ch15.ccd_diagram_v1 import CCD_energy, CCD_equations

for x in CCD_energy(2):
    print(f'E += {x}')

for x in CCD_equations(2):
    print(f'F2 += {x}')
