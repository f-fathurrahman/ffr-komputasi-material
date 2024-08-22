# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %% [markdown]
# # MintsHelper: Generating 1- and 2-electron Integrals with <span style='font-variant: small-caps'> Psi4 </span>
#
# In all of quantum chemistry, one process which is common to nearly every method is the evaluation of one-
# and two-electron integrals.  Fortunately, we can leverage infrastructure in <span style='font-variant: small-caps'> 
# Psi4 </span> to perform this task for us.  This tutorial will discuss the [``psi4.core.MintsHelper``](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper "Go to API") class, which is an
# interface for the powerful Psi4 ``libmints`` library which wraps the `libint` library, where these integrals are actually computed.  
#
# ## MintsHelper Overview
# In order to compute 1- and 2-electron integrals, we first need a molecule and basis set with which to work.  So, 
# before diving into `MintsHelper`, we need to build these objects.  In the cell below, we have imported
# <span style='font-variant: small-caps'> Psi4 </span> and NumPy, defined a water molecule, and set the basis to
# cc-pVDZ.  We've also set the memory available to <span style='font-variant: small-caps'> Psi4</span>, as well as
# defined a variable `numpy_memory` which we will discuss later.

# %%
# ==> Setup <==
# Import statements
import psi4
import numpy as np

# Memory & Output file
psi4.set_memory(int(2e9))
numpy_memory = 2
psi4.core.set_output_file('output.dat', False)

# Molecule definition
h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

# Basis Set
psi4.set_options({'basis': 'cc-pvdz'})

# %% [markdown]
# Now, we are ready to create an instance of the `MintsHelper` class.  To do this, we need to pass a `BasisSet`
# object to the `MintsHelper` initializer.  Fortunately, from the previous tutorial on the `Wavefunction` class, we know
# that we can obtain such an object from an existing wavefunction.  So, let's build a new wavefunction for our molecule,
# get the basis set object, and build an instance of `MintsHelper`:

# %%
# ==> Build MintsHelper Instance <==
# Build new wavefunction
wfn = psi4.core.Wavefunction.build(h2o, psi4.core.get_global_option('basis'))

# Initialize MintsHelper with wavefunction's basis set
mints = psi4.core.MintsHelper(wfn.basisset())

# %% [markdown]
# Below are summarized several commonly computed quantities and how to obtain them using a `MintsHelper` class method:
#
# | Quantity | Function | Description |
# |----------|----------|-------------|
# | AO Overlap integrals | [mints.ao_overlap()](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_overlap "Go to Documentation") | Returns AO overlap matrix as a `psi4.core.Matrix` object |
# | AO Kinetic Energy | [mints.ao_kinetic()](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_kinetic "Go to Documentation") | Returns AO kinetic energy matrix as a `psi4.core.Matrix` object |
# | AO Potential Energy | [mints.ao_potential()](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_potential "Go to Documentation") | Returns AO potential energy matrix as a `psi4.core.Matrix` object |
# | AO Electron Repulsion Integrals | [mints.ao_eri()](http://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_eri "Go to Documentation") | Returns AO electron repulsion integrals as a `psi4.core.Matrix` object 
#
# As discussed previously, any of these `psi4.core.Matrix` objects can be accessed as NumPy arrays, which is the preferred 
# method in Psi4NumPy.  For a Psi4 matrix `A`, we can access a NumPy view using `np.asarray(A)` or `A.np`, or we can make a
# copy of the matrix using `np.array(A)`.  This works as one would expect, converting square matrices into rank-2 NumPy 
# arrays, for the overlap (S), kinetic energy (T), and potential energy (V) matrices.  In Psi4, the electron repulsion integrals 
# (ERIs) are handled somewhat differently; `mints.ao_eri()` returns the rank-4 ERI tensor packed into a 2D matrix.  If the 
# four indices of the ERI are p, q, r, s, then this element of the Psi4 Matrix can be accessed by first computing composite 
# indices `pq = p * nbf + q` and `rs = r * nbf + s`, and then accessing element `A.get(pq,rs)`.  However, for convenience, 
# the NumPy view is a rank-4 tensor, and a particular ERI is more simply accessed like this:
# ~~~python
# I = np.asarray(mints.ao_eri())
# val = I[p][q][r][s]
# ~~~
#
# In addition to these methods, another which is worth mentioning is the [`MintsHelper.mo_eri()`](http://psicode.org
# /psi4manual/master/psi4api.html#psi4.core.MintsHelper.mo_eri "Go to Documentation") function, which can transform 
# the four-index, two-electron repulsion integrals from the atomic orbital (AO) to the molecular orbital (MO) basis,
# which will be important in MP2 theory.  
#
# ## Memory Considerations
#
# Before moving forward to computing any 1- or 2-electron integrals, we must first discuss the memory requirements of
# these objects.  Whenever these quantities are computed, they are stored directly in memory (a.k.a. RAM,
# *not* on the hard drive) which, for a typical laptop or personal computer, usually tops out at around 16 GB of 
# space.  The storage space required by the two-index AO overlap integrals and four-index ERIs scales as ${\cal O}(N^2)$ 
# and ${\cal O}(N^4)$, respectively, where $N$ is the number of AO basis functions.  This means that for a
# system with 500 AO basis functions, while the AO overlap integrals will only require 1 MB of memory to store,
# the ERIs will require a staggering **500 GB** of memory!! This can be reduced to **62.5 GB** of memory if integral permutational symmetry is used. 
# However, this complicates the bookkeeping, and is not used in the `mints` functions discussed above.  For this reason, as well as the steep computational 
# scaling of many of the methods demonstrated here, we limit ourselves to small systems ($\sim50$ basis functions)
# which should not require such egregious amounts of memory.  Additionally, we will employ a "memory check" to catch
# any case which could potentially try to use more memory than is available:
# ~~~python
# # Memory check for ERI tensor
# I_size = (nbf**4) * 8.e-9
# print('Size of the ERI tensor will be %4.2f GB.' % (I_size))
# memory_footprint = I_size * 1.5
# if I_size > numpy_memory:
#     psi4.core.clean()
#     raise Exception("Estimated memory utilization (%4.2f GB) exceeds allotted memory \
#                      limit of %4.2f GB." % (memory_footprint, numpy_memory))
# ~~~
# In this example, we have somewhat arbitrarily assumed that whatever other matrices we may need, in total their memory
# requirement will not exceed 50% of the size of the ERIs (hence, the total memory footprint of `I_size * 1.5`)
# Using the `numpy_memory` variable, we are able to control whether the ERIs will be computed, based on the amount of
# memory required to store them. 
#
# <font color="red">**NOTE: DO NOT EXCEED YOUR SYSTEM'S MEMORY.  THIS MAY RESULT IN YOUR PROGRAM AND/OR COMPUTER CRASHING!**</font>
#
# ## Examples: AO Overlap, AO ERIs, Core Hamiltonian
# The cell below demonstrates obtaining the AO overlap integrals, conducting the
# above memory check, and computing the ERIs and core Hamiltonian matrix for our water molecule.

# %%
# ==> Integrals galore! <==
# AO Overlap
S = np.asarray(mints.ao_overlap())

# Number of basis functions
nbf = S.shape[0]

# Memory check
I_size = (nbf ** 4) * 8.e-9
print('Size of the ERI tensor will be %4.2f GB.' % (I_size))
memory_footprint = I_size * 1.5
if I_size > numpy_memory:
    psi4.core.clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds allotted memory \
                     limit of %4.2f GB." % (memory_footprint, numpy_memory))

# Compute AO-basis ERIs
I = mints.ao_eri()

# Compute AO Core Hamiltonian
T = np.asarray(mints.ao_kinetic())
V = np.asarray(mints.ao_potential())
H = T + V
