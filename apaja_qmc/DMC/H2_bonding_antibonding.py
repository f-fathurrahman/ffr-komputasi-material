import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

zeta = 1.0
# proton positions
R_A = np.array([0.0,0.0,0.0])
R_B = np.array([1.4,0.0,0.0])

def phi_1s(rvec):
    r = np.sqrt((rvec**2).sum() )
    return exp(-zeta*r)

def sigma_g(r):
    r1 = r[0] 
    r2 = r[1]
    sigma_g = (phi_1s(r1-R_A) + phi_1s(r1-R_B)) * (phi_1s(r2-R_A) + phi_1s(r2-R_B))
    return sigma_g

def sigma_u(r):
    r1 = r[0] 
    r2 = r[1]
    sigma_u = (phi_1s(r1-R_A) - phi_1s(r1-R_B)) * (phi_1s(r2-R_A) - phi_1s(r2-R_B))
    return sigma_u


# Parameters
n_samples = 30000  # Number of samples
N = 2    # Number of electrons
dim = 3            # 3D system

e_positions = np.zeros((n_samples, N, dim))
r = np.random.random((N, dim))*5.0
rp = np.zeros((N, dim))
# *****************
# choose wave function
#psi = sigma_u
psi = sigma_g
# *****************
phi0 = psi(r)
step = 0.1
Ntry = 0
Nacc = 0

    


for sample in range(n_samples):
    print('sample',sample)
    if sample==0:
        rr = 10000 # thermalize
    else:
        rr = 3
        
    for i in range(rr):
        d = np.random.normal(size=(N,dim))*step
        rp[:] = r+d
        phi = psi(rp)
        Ntry+=1
        accept = False
        ratio = phi**2/phi0**2
        if ratio>1.0:
            accept=True
        elif ratio>np.random.random():
            accept =  True

        if accept:
            phi0 = phi
            r[:] = rp
            Nacc += 1

        if Ntry%100:
            acceptance = Nacc*100.0/Ntry
            if acceptance>60.0:
                step *= 1.1
            elif  acceptance<50.0:
                step *= 0.9

    e_positions[sample] = r

# save samples to file
# save coordinates as x y z particle_index (all as floats)  
p_indices = np.repeat(np.arange(N), n_samples) # repeat 0.0 and 1.0
resh = e_positions.reshape(-1, 3)
comb_ps = np.column_stack((resh, p_indices))
with open('tt','w') as f:    
    np.savetxt(f, comb_ps, fmt='%.6f')

print('data in file tt')
print('plot with gnuplot, type in gnuplot prompt:')
print('unset key')
print("sp [-3:3][-3:3][-3:3] 'tt' u 1:2:3:((int($4)==0?6:7))  w p ps 0.001  lc variable")
# protons:
# rep "< echo '1.4, 0, 0'" w p ps 2 pt 7 lc 8
# rep "< echo '0, 0, 0'" w p ps 2 pt 7 lc 8



