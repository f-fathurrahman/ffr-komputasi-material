"""
compare_fullgrid_vs_ci_3e_1d_debug.py
------------------------------------
Debug version for N=10, M=4, 3-electron 1D system.

Purpose:
- Verify consistency between full-grid direct discretization and CI energies.
- Explicitly test normalization, antisymmetry, and reconstruction.

Run: python compare_fullgrid_vs_ci_3e_1d_debug.py
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import itertools, math

# ---------------- Parameters ----------------
N = 20
x_min, x_max = -8, 8
x = np.linspace(x_min, x_max, N)
dx = x[1] - x[0]
dx3 = dx**3
m = 1.0
hbar = 1.0
omega = 1.0
a_soft = 1.0
M = 5   # number of spatial orbitals

# ---------------- Single-particle Hamiltonian ----------------
diag = -2.0*np.ones(N)/dx**2
off  = 1.0*np.ones(N-1)/dx**2
L1 = sp.diags([off, diag, off], [-1,0,1], format='csr')
Vx = 0.5*m*(omega**2)*x**2
H1 = - (hbar**2)/(2*m)*L1 + sp.diags(Vx,0)

eps, phi = spla.eigsh(H1, k=M, which='SA')
eps, phi = np.real(eps), np.real(phi)
idx = np.argsort(eps)
eps, phi = eps[idx], phi[:,idx]
for i in range(M):
    phi[:,i] /= math.sqrt(np.sum(phi[:,i]**2)*dx)

print("Single-particle energies:", np.round(eps,6))

# ---------------- Two-electron integrals ----------------
X = x.reshape(N,1)
Vpair = 1.0/np.sqrt((X-X.T)**2 + a_soft**2)
h = np.zeros((M,M))
g = np.zeros((M,M,M,M))
for p in range(M):
    for q in range(M):
        h[p,q] = np.sum(phi[:,p]*(H1@phi[:,q]))*dx
for p in range(M):
    for q in range(M):
        for r in range(M):
            for s in range(M):
                integrand = (phi[:,p]*phi[:,r])[:,None]*Vpair*(phi[:,q]*phi[:,s])[None,:]
                g[p,q,r,s] = np.sum(integrand)*dx**2

print("h matrix:\n", np.round(h,4))

# ---------------- Build CI (all-spin-up) ----------------
from itertools import combinations
dets = list(combinations(range(M),3))
# ---------- Replace CI building with Slater-Condon based builder ----------
# dets is list of tuples of 3 spatial orbital indices (sorted)
D = len(dets)
Hci = np.zeros((D, D))

def insertion_parity_remove_insert(base_det, removed, added):
    """
    Compute phase = (-1)^{#swaps} for going from base_det to new_det where:
    - base_det: tuple/list (sorted) representing determinant (J)
    - removed: single index removed from base_det (value)
    - added: single index to add (value)
    Returns parity_sign (±1).
    """
    cur = list(base_det)
    # index of removed element in current ordering
    idx_remove = cur.index(removed)
    cur.pop(idx_remove)  # list with one element removed
    # find insertion position for 'added' to keep ordering
    pos = 0
    while pos < len(cur) and cur[pos] < added:
        pos += 1
    # total swaps = idx_remove + pos
    return (-1) ** (idx_remove + pos)

def double_exc_parity(base_det, removed_pair, added_pair):
    """
    Compute parity for double excitation: remove two indices (in the order they
    appear in base_det) then insert two new indices so the new list is sorted.
    Returns ±1.
    """
    cur = list(base_det)
    # remove in decreasing index order of the positions so pop doesn't shift earlier positions wrongly
    # but the parity depends on original positions; do as follows:
    # find indices (positions) of the two removed orbitals in base_det
    pos_r = [cur.index(removed_pair[0]), cur.index(removed_pair[1])]
    # sort by descending position to pop without reindexing earlier one wrongly
    removed_sorted_by_pos = [rp for _, rp in sorted(zip(pos_r, removed_pair), reverse=True)]
    total_swaps = 0
    for r in removed_sorted_by_pos:
        pos_rm = cur.index(r)
        total_swaps += pos_rm
        cur.pop(pos_rm)
    # now cur has length 1, insert the two added orbitals in sorted order, count insertions
    # for each insertion, count how many current entries are less than new added (that's insertion pos)
    added_sorted = sorted(list(added_pair))
    for a in added_sorted:
        pos = 0
        while pos < len(cur) and cur[pos] < a:
            pos += 1
        total_swaps += pos
        cur.insert(pos, a)
    return (-1) ** total_swaps

# Build Hci using Slater-Condon rules for spatial determinants
for I_idx, I in enumerate(dets):
    setI = set(I)
    for J_idx, J in enumerate(dets):
        setJ = set(J)
        # difference sets
        removed = sorted(list(setJ - setI))
        added = sorted(list(setI - setJ))
        ndiff = len(removed)  # must equal len(added)
        if ndiff == 0:
            # diagonal: sum of one-body plus 1/2 sum_{ij in occ} (g_{ijij} - g_{ijji})
            E1 = sum(h[p,p] for p in I)
            E2 = 0.0
            for p in I:
                for q in I:
                    E2 += 0.5 * (g[p,q,p,q] - g[p,q,q,p])
            Hci[I_idx, J_idx] = E1 + E2
        elif ndiff == 1:
            # single excitation: I = (J with q->p)
            p = added[0]
            q = removed[0]
            # parity sign from removing q and inserting p into J
            phase = insertion_parity_remove_insert(J, q, p)
            # one-body term + sum over occupied r in intersection
            common = sorted(list(setI.intersection(setJ)))
            val = h[p, q]
            for r in common:
                val += (g[p, r, q, r] - g[p, r, r, q])
            Hci[I_idx, J_idx] = phase * val
        elif ndiff == 2:
            # double excitation: I = (J with r,s -> p,q)
            p, q = added
            r, s = removed
            # parity sign for double excitation
            phase = double_exc_parity(J, removed, added)
            # matrix element = (g_{p q r s} - g_{p q s r})
            Hci[I_idx, J_idx] = phase * (g[p, q, r, s] - g[p, q, s, r])
        else:
            # triple or higher difference -> zero
            Hci[I_idx, J_idx] = 0.0





Eci, Cci = np.linalg.eigh(Hci)
print("\nCI energies (M=4):", np.round(Eci,6))
ci0 = Cci[:,0]
print("CI coefficients:", np.round(ci0,3))

# ---------------- Full-grid Hamiltonian ----------------
I = sp.eye(N, format='csr')
T3 = - (hbar**2)/(2*m) * (sp.kron(sp.kron(L1,I),I) + sp.kron(sp.kron(I,L1),I) + sp.kron(sp.kron(I,I),L1))
Vext = (np.kron(np.kron(Vx, np.ones(N)), np.ones(N))
       + np.kron(np.kron(np.ones(N),Vx), np.ones(N))
       + np.kron(np.kron(np.ones(N),np.ones(N)),Vx))

# interactions on diagonal
Vint = np.zeros((N,N,N))
for i in range(N):
    for j in range(N):
        for k in range(N):
            Vint[i,j,k] = (Vpair[i,j]+Vpair[i,k]+Vpair[j,k])
Vfull = Vext + Vint.ravel()
Hfull = T3 + sp.diags(Vfull,0)

print(f"\nFull-grid Hamiltonian size: {N**3}")

vals, vecs = spla.eigsh(Hfull, k=6, which='SA')
print("Raw full-grid eigenvalues:", np.round(vals,6))

# ---------------- Antisymmetrize ----------------
def antisymmetrize(P):
    return (P - P.transpose(1,0,2) - P.transpose(2,1,0)
              - P.transpose(0,2,1) + P.transpose(1,2,0) + P.transpose(2,0,1))/6.0

def rayleigh(H,psi):
    return (np.vdot(psi,H@psi)*dx3/np.vdot(psi,psi)/dx3).real

bestE = 1e9
bestpsi = None
print("len(vals) = ", len(vals))
for i in range(len(vals)):
    psi = vecs[:,i].reshape(N,N,N)
    psiA = psi #antisymmetrize(psi)
    normA = np.sqrt(np.vdot(psiA.ravel(),psiA.ravel())*dx3)
    print("normA = ", normA)
    if normA < 1e-10:
        continue
    psiA = psiA.ravel()/normA
    E = rayleigh(Hfull,psiA)
    print(f"Eigen {i}: eigval={vals[i]:.6f}, antisymm E={E:.6f}")
    if E < bestE:
        bestE, bestpsi = E, psiA

print(f"\nBest antisymmetric energy: {bestE:.6f}")

# ---------------- Compare CI and full-grid ----------------
print(f"CI ground energy = {Eci[0]:.6f},  Full-grid antisymm = {bestE:.6f},  Δ = {bestE-Eci[0]:.3e}")

exit()


# ---------- Diagnostic 1: Inspect direct/exchange integrals ----------
print("\\n--- g direct/exchange checks (sample pairs) ---")
for p in range(min(4,M)):
    for q in range(min(4,M)):
        direct = g[p,q,p,q]   # (p q|p q)  often written g_{pqpq}
        exchange = g[p,q,q,p] # (p q|q p)  often written g_{pqqp}
        print(f"p={p}, q={q}: direct g[{p},{q},{p},{q}] = {direct:.6e}, exchange g[{p},{q},{q},{p}] = {exchange:.6e}")
print("Note: direct integrals (p==q or p!=q) should be >= 0 (Coulomb is repulsive).")

# ---------- Diagnostic 2: Compare Hci diagonal vs Slater-Condon explicit formula ----------
print("\\n--- Hci diagonal vs Slater-Condon formula ---")
Hci_diag = np.diag(Hci).copy()
for I_idx, det in enumerate(dets):
    # det contains spatial orbital indices (0..M-1)
    # one-body sum
    E1 = sum(h[p,p] for p in det)
    # two-body sum (i,j over occupied)
    E2 = 0.0
    for i in det:
        for j in det:
            E2 += 0.5 * (g[i,j,i,j] - g[i,j,j,i])
    E_SC = E1 + E2
    print(f"Det {I_idx} {det}: Hci_diag={Hci_diag[I_idx]:.8f}, Slater-Condon={E_SC:.8f} (E1={E1:.6f}, E2={E2:.6f})")
# If these numbers differ, the CI-builder (apply_ops loops) is wrong.

# ---------- Diagnostic 3: Print full H expectation for single-determinant grid reconstruction ----------
# Reconstruct spatial wavefunction for first determinant det0 = dets[0] (orbitals p,q,r).
print("\\n--- Reconstruct single-determinant |p q r> on grid and evaluate <H> ---")
det0 = dets[0]
print("det0 spatial orbitals:", det0)
psi_det0 = np.zeros((N,N,N), dtype=float)
# build determinant value at each grid point (i,j,k)
for i in range(N):
    for j in range(N):
        for k in range(N):
            Mmat = np.array([
                [phi[i, det0[0]], phi[i, det0[1]], phi[i, det0[2]]],
                [phi[j, det0[0]], phi[j, det0[1]], phi[j, det0[2]]],
                [phi[k, det0[0]], phi[k, det0[1]], phi[k, det0[2]]]
            ], dtype=float)
            psi_det0[i,j,k] = np.linalg.det(Mmat)
psi_det0_flat = psi_det0.ravel()
norm_det0 = math.sqrt(np.vdot(psi_det0_flat, psi_det0_flat) * dx3)
print("norm(det0) (with dx^3) =", norm_det0)
psi_det0_flat /= norm_det0
E_det0_via_H = (np.vdot(psi_det0_flat, (Hfull.dot(psi_det0_flat))) * dx3 / (np.vdot(psi_det0_flat, psi_det0_flat) * dx3)).real
# Compare to Slater-Condon expectation for the same determinant (E_SC computed above)
# find index 0 in previous loop: ensure E_SC_for_det0 is available
E1 = sum(h[p,p] for p in det0)
E2 = 0.0
for i in det0:
    for j in det0:
        E2 += 0.5 * (g[i,j,i,j] - g[i,j,j,i])
E_SC_det0 = E1 + E2
print(f"E_det0 (grid eval) = {E_det0_via_H:.8f}, Slater-Condon = {E_SC_det0:.8f}, (E1={E1:.6f}, E2={E2:.6f})")

# ---------- Diagnostic 4: Reconstruct entire CI ground state on grid and evaluate <H> ----------
print("\\n--- Reconstruct CI ground state on grid and evaluate <H> ---")
ci0 = Cci[:,0]   # CI ground coefficients
psi_ci = np.zeros((N,N,N), dtype=float)
for coeff, det in zip(ci0, dets):
    # det is tuple of 3 spatial orbital indices
    for i in range(N):
        for j in range(N):
            for k in range(N):
                Mmat = np.array([
                    [phi[i, det[0]], phi[i, det[1]], phi[i, det[2]]],
                    [phi[j, det[0]], phi[j, det[1]], phi[j, det[2]]],
                    [phi[k, det[0]], phi[k, det[1]], phi[k, det[2]]]
                ], dtype=float)
                psi_ci[i,j,k] += coeff * np.linalg.det(Mmat)
psi_ci_flat = psi_ci.ravel()
norm_ci = math.sqrt(np.vdot(psi_ci_flat, psi_ci_flat) * dx3)
psi_ci_flat /= norm_ci
# energy via Hfull
E_ci_via_H = (np.vdot(psi_ci_flat, (Hfull.dot(psi_ci_flat))) * dx3 / (np.vdot(psi_ci_flat, psi_ci_flat) * dx3)).real
print(f"CI eigval (from diagonalization) = {Eci[0]:.8f}")
print(f"CI reconstructed energy via full H = {E_ci_via_H:.8f}")

# Overlap between CI reconstructed and first antisymmetric full-grid state (if you saved it)
# If 'bestpsi' available from earlier, compute:
#try:
#    ov = abs(np.vdot(bestpsi, psi_ci_flat)) * dx3
#    print(f"Overlap with best full-grid antisymm = {ov:.8e}")
#except NameError:
#    pass
