# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import ase.io
from mini_mlip import DescriptorSO3

# %%
import mini_mlip

# %%
import numpy as np

# %%
from scipy.special import sph_harm_y, spherical_in

# %%
atoms_list = ase.io.read("DATASET_OTHERS/TiAl_gabung.xyz@:")

# %%
atoms = atoms_list[0]

# %%
atoms.get_potential_energy()

# %%
atoms.get_forces()

# %%
lmax = 4
nmax = 3
rcut = 3.5
α = 2.0
SO3_calc = DescriptorSO3(
    nmax=nmax,
    lmax=lmax,
    rcut=rcut,
    alpha=α,
    derivative=True,
    stress=False,
    cutoff_function='cosine'
)

# %%
# Set a copy ?
SO3_calc._atoms = atoms # ._atoms is not referenced anymore

# %%
SO3_calc.build_neighbor_list() # atoms_ids is None

# %% [markdown]
# The following variables are initialized after call to build_neighbor_list:

# %%
SO3_calc.center_atoms.shape

# %%
SO3_calc.neighborlist.shape

# %%
SO3_calc.seq

# %%
SO3_calc.atomic_weights.shape

# %%
SO3_calc.neighbor_indices.shape

# %% [markdown]
# This will initialize more arrays:

# %%
SO3_calc.initialize_arrays() # intialize _plist, _dplist, _pstress

# %%
SO3_calc._plist.shape

# %%
SO3_calc._dplist.shape

# %%
SO3_calc._pstress.shape

# %%
ncoefs = SO3_calc.nmax*(SO3_calc.nmax+1)//2*(SO3_calc.lmax+1)

# %%
ncoefs

# %%
tril_indices = np.tril_indices(SO3_calc.nmax, k=0)

# %%
tril_indices[0]

# %%
tril_indices[1]

# %%
ls = np.arange(SO3_calc.lmax+1)

# %%
# ls

# %%
norm = np.sqrt(2*np.sqrt(2)*np.pi/np.sqrt(2*ls+1))

# %%
norm


# %% [markdown]
# This call is quite complicated:

# %%
def my_compute_cs(pos, nmax, lmax, rcut, alpha, cutoff):

    # compute the overlap matrix
    w = mini_mlip.W(nmax)

    # get the norm of the position vectors
    Ris = np.linalg.norm(pos, axis=1) # (Nneighbors)

    # initialize Gauss Chebyshev Quadrature
    GCQuadrature, weight = mini_mlip.GaussChebyshevQuadrature(nmax,lmax) #(Nquad)
    weight *= rcut/2
    # transform the quadrature from (-1,1) to (0, rcut)
    Quadrature = rcut/2*(GCQuadrature+1)

    # compute the arguments for the bessel functions
    BesselArgs = 2*alpha*np.outer(Ris,Quadrature)#(Nneighbors x Nquad)

    # initalize the arrays for the bessel function values
    # and the G function values
    Bessels = np.zeros((len(Ris), len(Quadrature), lmax+1), dtype=np.float64) #(Nneighbors x Nquad x lmax+1)
    Gs = np.zeros((nmax, len(Quadrature)), dtype=np.float64) # (nmax, nquad)

    # compute the g values
    for n in range(1,nmax+1,1):
        Gs[n-1,:] = mini_mlip.g(Quadrature, n, nmax, rcut, w)

    # compute the bessel values
    for l in range(lmax+1):
        Bessels[:,:,l] = spherical_in(l, BesselArgs)

    # mutliply the terms in the integral separate from the Bessels
    Quad_Squared = Quadrature**2
    Gs *= Quad_Squared * np.exp(-alpha*Quad_Squared) * np.sqrt(1-GCQuadrature**2) * weight

    # perform the integration with the Bessels
    integral_array = np.einsum('ij,kjl->kil', Gs, Bessels) # (Nneighbors x nmax x lmax+1)

    # compute the gaussian for each atom and multiply with 4*pi
    # to minimize floating point operations
    # weight can also go here since the Chebyshev gauss quadrature weights are uniform
    exparray = 4*np.pi*np.exp(-alpha*Ris**2) # (Nneighbors)

    cutoff_array = cutoff(Ris, rcut) # cutoff is a Function

    exparray *= cutoff_array

    # get the spherical coordinates of each atom
    thetas = np.arccos(pos[:,2]/Ris[:])
    phis = np.arctan2(pos[:,1], pos[:,0])

    # determine the size of the m axis
    msize = 2*lmax+1
    # initialize an array for the spherical harmonics
    ylms = np.zeros((len(Ris), lmax+1, msize), dtype=np.complex128)

    # compute the spherical harmonics
    for l in range(lmax+1):
        for m in range(-l,l+1,1):
            midx = msize//2 + m
            ylms[:,l,midx] = sph_harm_y(m, l, phis, thetas)

    # multiply the spherical harmonics and the radial inner product
    Y_mul_innerprod = np.einsum('ijk,ilj->iljk', ylms, integral_array)

    # multiply the gaussians into the expression
    C = np.einsum('i,ijkl->ijkl', exparray, Y_mul_innerprod)
    return C


# %%
SO3_calc.neighborlist.shape

# %%
SO3_calc.primitive

# %%
SO3_calc.rcut

# %%
SO3_calc.neighborlist

# %%
type(SO3_calc.nmax)

# %%
cs = my_compute_cs(
    SO3_calc.neighborlist,
    SO3_calc.nmax,
    SO3_calc.lmax,
    SO3_calc.rcut,
    SO3_calc.alpha,
    SO3_calc._cutoff_function
)
cs *= SO3_calc.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis]

# %%
cs.shape

# %%
cs.shape

# %%
cs = np.einsum('inlm,l->inlm', cs, norm)
# everything good up to here
for i in np.unique(SO3_calc.seq[:,0]):
    centers = SO3_calc.neighbor_indices[:,0] == i
    ctot = cs[centers].sum(axis=0)
    P = np.einsum('ijk,ljk->ilj', ctot, np.conj(ctot)).real
    SO3_calc._plist[i] = P[tril_indices].flatten()

# %%
desc_vector_dict = {
    'x' : SO3_calc._plist,
    'dxdr' : None,
     'rdxdr' : None,
     'elements' : list(atoms.symbols)
}

# %%
desc_vector_dict["x"].shape

# %%
