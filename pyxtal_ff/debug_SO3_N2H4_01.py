import ase.io
from mini_mlip import DescriptorSO3, compute_dcs, compute_cs
import numpy as np


def calculate_SO3_power_spectrum(SO3_calc, atoms, atom_ids=None):

    # Set a copy ?
    SO3_calc._atoms = atoms # ._atoms is not referenced anymore ?

    SO3_calc.build_neighbor_list(atom_ids)
    SO3_calc.initialize_arrays()

    ncoefs = SO3_calc.nmax*(SO3_calc.nmax+1)//2*(SO3_calc.lmax+1)
    tril_indices = np.tril_indices(SO3_calc.nmax, k=0)

    ls = np.arange(SO3_calc.lmax+1)
    norm = np.sqrt(2*np.sqrt(2)*np.pi/np.sqrt(2*ls+1))

    if SO3_calc.derivative:
        # get expansion coefficients and derivatives
        cs, dcs = compute_dcs(
            SO3_calc.neighborlist,
            SO3_calc.nmax,
            SO3_calc.lmax,
            SO3_calc.rcut,
            SO3_calc.alpha,
            SO3_calc._cutoff_function
        )
        # weight cs and dcs
        cs *= SO3_calc.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis]
        dcs *= SO3_calc.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
        cs = np.einsum('inlm,l->inlm', cs, norm)
        dcs = np.einsum('inlmj,l->inlmj', dcs, norm)
        Ris = SO3_calc.center_atoms
        Rjs = SO3_calc.neighborlist + Ris
        for i in np.unique(SO3_calc.seq[:,0]):
            # find atoms for which i is the center
            centers = SO3_calc.neighbor_indices[:,0] == i
            # find neighbors for which i is not the index
            neighs = SO3_calc.neighbor_indices[:,1] != i
            # get the indices for both conditions
            inds = centers*neighs
            # total up the c array for the center atom
            ctot = cs[centers].sum(axis=0)
            #ctot = np.einsum('nlm,l->nlm', ctot,norm)
            # get dc weights
            # compute the power spectrum
            P = np.einsum('ijk,ljk->ilj', ctot, np.conj(ctot)).real
            # compute the gradient of the power spectrum for each neighbor
            dP = np.einsum('wijkn,ljk->wiljn', dcs[centers], np.conj(ctot))
            dP += np.conj(np.transpose(dP, axes=[0,2,1,3,4]))
            dP = dP.real

            rdPi = np.einsum('wn,wijkm->wijknm', Ris[centers], dP)
            rdPj = np.einsum('wn,wijkm->wijknm', Rjs[centers], dP)
            # get ij pairs for center atom
            ijs = SO3_calc.neighbor_indices[centers]
            # loop over unique neighbor indices
            for j in np.unique(ijs[:,1]):
                # get the location of ij pairs in the NL
                # and therefore dP
                ijlocs = SO3_calc.neighbor_indices[centers,1] == j
                # get the location of the dplist element
                temp = SO3_calc.seq == np.array([i,j])
                seqloc = temp[:,0]*temp[:,1]
                # sum over ij pairs
                dPsum = np.sum(dP[ijlocs], axis=0)
                rdPjsum = np.sum(rdPj[ijlocs], axis=0)
                # flatten into dplist and rdplist
                SO3_calc._dplist[seqloc] += (dPsum[tril_indices].flatten()).reshape(ncoefs,3)
                SO3_calc._pstress[seqloc] -= (rdPjsum[tril_indices].flatten()).reshape(ncoefs,3,3)

            # get unique elements and store in feature vector
            SO3_calc._plist[i] = P[tril_indices].flatten()
            # get location if ii pair in seq
            temp = SO3_calc.seq == np.array([i,i])
            iiloc = temp[:,0]*temp[:,1]
            # get location of all ijs in seq
            ilocs = SO3_calc.seq[:,0] == i
            SO3_calc._dplist[iiloc] -= np.sum(SO3_calc._dplist[ilocs],axis=0)
            rdPisum = np.sum(rdPi, axis=0)
            SO3_calc._pstress[iiloc] += (rdPisum[tril_indices].flatten()).reshape(ncoefs,3,3)


        x = {
            'x' : SO3_calc._plist,
            'dxdr' : SO3_calc._dplist,
            'elements' : list(atoms.symbols),
            'seq' : SO3_calc.seq
            }
        if SO3_calc._stress:
            vol = atoms.get_volume()
            x['rdxdr'] = -SO3_calc._pstress/vol
        else:
            x['rdxdr'] = None

    else:
        cs = compute_cs(
            SO3_calc.neighborlist,
            SO3_calc.nmax,
            SO3_calc.lmax,
            SO3_calc.rcut,
            SO3_calc.alpha,
            SO3_calc._cutoff_function
        )
        cs *= SO3_calc.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis]
        cs = np.einsum('inlm,l->inlm', cs, norm)
        # everything good up to here
        for i in np.unique(SO3_calc.seq[:,0]):
            centers = SO3_calc.neighbor_indices[:,0] == i
            ctot = cs[centers].sum(axis=0)
            P = np.einsum('ijk,ljk->ilj', ctot, np.conj(ctot)).real
            SO3_calc._plist[i] = P[tril_indices].flatten()
        x = {'x' : SO3_calc._plist,
             'dxdr' : None,
             'rdxdr' : None,
             'elements' : list(atoms.symbols)
            }

    SO3_calc.clear_memory()
    return x


#atoms = ase.io.read("DATASET_N2H4_v2/N2H4_2mol_1data.xyz")


atoms_list = ase.io.read("DATASET_OTHERS/TiAl_gabung.xyz@:")
atoms = atoms_list[0]


# Should be invariant with w.r.t translations
#pos_shifted = atoms.positions.copy()
#pos_shifted[:,2] = atoms.positions[:,2] + 3.0
#atoms.set_positions(pos_shifted)

lmax = 4
nmax = 3
rcut = 3.5
alpha = 2.0

desc_calc = DescriptorSO3(
    nmax=nmax,
    lmax=lmax,
    rcut=rcut,
    alpha=alpha,
    derivative=True,
    stress=False,
    cutoff_function='cosine'
)
x = desc_calc.calculate(atoms)

print(desc_calc)
print(x["x"][0])
print(x["x"][1])
print(x["x"].shape)