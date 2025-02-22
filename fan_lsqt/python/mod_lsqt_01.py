import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy import special

def find_index(nx, ny, Ny):
    index = nx * Ny + ny
    return index

def find_H(Nx, Ny, W):
    N = Nx * Ny; # total number of sites
    row_H = np.zeros(N * 4 - Nx * 2) # the row indices for H
    col_H = np.zeros(N * 4 - Nx * 2) # the column indices for H
    Hij = -np.ones(N * 4 - Nx * 2, complex) # nonzero Hamiltonian elements
    row_V = np.zeros(N * 2) # row indices for V
    col_V = np.zeros(N * 2) # column indices for V
    Vij = np.zeros(N * 2, complex) # nonzero velocity matrix elements
    row_U = np.arange(0, N) # row and column indices for U
    Uij = np.random.uniform(-W * 0.5, W * 0.5, N) # on-site potentials
    count_H = 0 # number of nonzero H elements
    count_V = 0 # number of nonzero V elements
    for nx in range(Nx):
        for ny in range(Ny):
            # (0) # get the index of the center site
            index_center = find_index(nx, ny, Ny)
            # (1) consider the left neighbor (periodic boundary)
            index_left = find_index((nx - 1) % Nx, ny, Ny)
            row_H[count_H] = index_center
            col_H[count_H] = index_left
            count_H += 1
            row_V[count_V] = index_center
            col_V[count_V] = index_left
            Vij[count_V] = 1j
            count_V += 1
            # (2) consider the right neighbor (periodic boundary)
            index_right = find_index((nx + 1) % Nx, ny, Ny)
            row_H[count_H] = index_center
            col_H[count_H] = index_right
            count_H += 1
            row_V[count_V] = index_center
            col_V[count_V] = index_right
            Vij[count_V] = -1j
            count_V += 1
            # (3) consider the upper neighbor (open boundary)
            if ny < Ny - 1:
                index_up = find_index(nx, (ny + 1), Ny)
                row_H[count_H] = index_center
                col_H[count_H] = index_up
                count_H += 1
            # (4) consider the down neighbor (open boundary)
            if ny > 0:
                index_down = find_index(nx, (ny - 1), Ny)
                row_H[count_H] = index_center
                col_H[count_H] = index_down
                count_H += 1
    H = sparse.csr_matrix((Hij, (row_H, col_H)), shape = (N, N))
    U = sparse.csr_matrix((Uij, (row_U, row_U)), shape = (N, N))
    H = H + U
    V = sparse.csr_matrix((Vij, (row_V, col_V)), shape = (N, N))
    return (H, V)

def create_state(N):
    random_phase = np.random.uniform(0, 2 * np.pi, N)
    phi = np.cos(random_phase) + np.sin(random_phase) * 1j; 
    phi /= np.linalg.norm(phi);
    return phi


def find_moments(M, H_scaled, phi_left, phi_right):  
    C = np.zeros(M)
    phi_0 = phi_right
    C[0] = np.vdot(phi_left, phi_0).real
    phi_1 = H_scaled.dot(phi_0)
    C[1] = np.vdot(phi_left, phi_1).real
    for m in range(2, M):
        phi_2 = 2.0 * H_scaled.dot(phi_1) - phi_0
        C[m] = np.vdot(phi_left, phi_2).real
        phi_0 = phi_1
        phi_1 = phi_2
    return C


def jackson_damping(M):
    m = np.arange(M)
    a = 1.0 / (M + 1.0)
    g = (1.0 - m * a) * np.cos(np.pi * m * a)
    g += a * np.sin(np.pi * m * a) / np.tan(np.pi * a)
    return g


def plot_jackson_damping():
    M = 10000
    plt.plot(np.arange(M), jackson_damping(M))
    plt.xlabel('m')
    plt.ylabel('g')
    plt.show()


def chebyshev_summation(M, C, E_scaled, E_max):
    g = jackson_damping(M) # get the damping factor
    C *= g                 # apply the damping factor
    Ne = E_scaled.size     # number of energy points
    T0 = np.ones(Ne)
    T1 = E_scaled
    s = C[1] * T1
    for m in range(2, M):
        T2 = 2.0 * E_scaled * T1 - T0
        T0 = T1
        T1 = T2
        s += C[m] * T2
    s *= 2.0
    s += C[0] * np.ones(Ne)
    s *= 2.0 / (np.pi * np.sqrt(1.0 - E_scaled * E_scaled) * E_max)
    return s


def find_dos(M, E_max, E_scaled, H_scaled, phi):
    C = find_moments(M, H_scaled, phi, phi)
    DOS = chebyshev_summation(M, C, E_scaled, E_max)
    return DOS

def evolve(H_scaled, dt_scaled, sign, phi_i):
    phi_0 = phi_i
    phi_1 = H_scaled.dot(phi_i)
    j0 = special.jv(0, dt_scaled)
    j1 = special.jv(1, dt_scaled)
    phi_o = j0 * phi_0 + 2.0 * (-1j * sign) * j1 * phi_1
    m = 2
    while True:
        jm = special.jv(m, dt_scaled)
        if abs(jm) < 1.0e-15:
            break
        phi_2 = 2.0 * H_scaled.dot(phi_1) - phi_0
        phi_o += 2.0 * (-1j * sign) ** m * jm * phi_2
        phi_0 = phi_1
        phi_1 = phi_2
        m += 1
    return phi_o


def evolvex(H_scaled, V_scaled, dt_scaled, phi_i):
    phi_0 = phi_i
    phix_0 = phi_i * 0.0
    phi_1 = H_scaled.dot(phi_0)
    phix_1 = 1j * V_scaled.dot(phi_0)
    phi_o = (-2.0j) * special.jv(1, dt_scaled) * phix_1
    m = 2
    while True:
        jm = special.jv(m, dt_scaled)
        if abs(jm) < 1.0e-15:
            break
        phi_2 = 2.0 * H_scaled.dot(phi_1) - phi_0
        phix_2 = 2.0j * V_scaled.dot(phi_1)
        phix_2 += 2.0 * H_scaled.dot(phix_1) - phix_0
        phi_o += 2.0 * (-1j) ** m * jm * phix_2
        phi_0 = phi_1 
        phi_1 = phi_2
        phix_0 = phix_1
        phix_1 = phix_2
        m += 1
    return phi_o


def find_vac(M, E_max, dt_scaled, E_scaled, H_scaled, V, phi, DOS):
    Ne = E_scaled.size  # number of energy points
    Nt = dt_scaled.size # number of time steps
    dt = dt_scaled / E_max
    phi_left = phi
    phi_right = V.dot(phi)
    DOS_times_VAC_old = np.zeros(Ne)
    DOS_times_VAC_new = np.zeros(Ne)
    VAC = np.zeros((Nt, Ne))
    sigma_from_VAC = np.zeros((Nt, Ne))
    for nt in range(Nt):
        C = find_moments(M, H_scaled, V.dot(phi_left), phi_right)
        DOS_times_VAC_new = chebyshev_summation(M, C, E_scaled, E_max)
        VAC[nt, :] = DOS_times_VAC_new / DOS
        if nt > 0:
            tmp = dt[nt - 1] * (DOS_times_VAC_old + DOS_times_VAC_new) * 0.5
            sigma_from_VAC[nt, :] = sigma_from_VAC[nt - 1, :] + tmp
        DOS_times_VAC_old = DOS_times_VAC_new
        phi_left = evolve(H_scaled, dt_scaled[nt], -1, phi_left)
        phi_right = evolve(H_scaled, dt_scaled[nt], -1, phi_right)
    sigma_from_VAC *= 2.0 * np.pi # from e^2/hbar to e^2/h
    return VAC, sigma_from_VAC



def find_msd(M, E_max, dt_scaled, E_scaled, H_scaled, V_scaled, phi, DOS):
    Ne = E_scaled.size  # number of energy points
    Nt = dt_scaled.size # number of time steps
    dt = dt_scaled / E_max
    MSD = np.zeros((Nt, Ne))
    sigma_from_MSD = np.zeros((Nt, Ne))
    phix = phi * 0.0
    DOS_times_MSD_old = np.zeros(Ne)
    DOS_times_MSD_new = np.zeros(Ne)
    for nt in range(Nt):
        phix = evolve(H_scaled, dt_scaled[nt], 1, phix);
        phix += evolvex(H_scaled, V_scaled, dt_scaled[nt], phi);
        phi = evolve(H_scaled, dt_scaled[nt], 1, phi);
        C = find_moments(M, H_scaled, phix, phix);
        DOS_times_MSD_new = chebyshev_summation(M, C, E_scaled, E_max);
        MSD[nt, :] = DOS_times_MSD_new / DOS
        sigma_from_MSD[nt, :] = (DOS_times_MSD_new - DOS_times_MSD_old) / dt[nt] * 0.5
        DOS_times_MSD_old = DOS_times_MSD_new
    sigma_from_MSD *= 2.0 * np.pi # from e^2/hbar to e^2/h
    return MSD, sigma_from_MSD



def lsqt(Nx, Ny, W, M, E_max, E, dt):
    
    print("Find H ...", end=" ")
    H, V = find_H(Nx, Ny, W)
    print("... end of find H")
    
    phi = create_state(Nx * Ny)
    
    print("Find DOS ...", end=" ")
    DOS = find_dos(M, E_max, E/E_max, H/E_max, phi)
    print("... end of find DOS")

    print("Find VAC ... ", end=" ")
    VAC, sigma_from_VAC = find_vac(M, E_max, dt*E_max, E/E_max, H/E_max, V, phi, DOS)
    print("... end of find VAC")

    print("Finding MSD ... ", end=" ")
    MSD, sigma_from_MSD = find_msd(M, E_max, dt*E_max, E/E_max, H/E_max, V/E_max, phi, DOS)
    print("... end of find MSD")

    return DOS, VAC, sigma_from_VAC, MSD, sigma_from_MSD