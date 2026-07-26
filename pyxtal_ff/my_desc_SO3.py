from ase.neighborlist import NeighborList, PrimitiveNeighborList
import numpy as np
from scipy.special import sph_harm_y, spherical_in

def Cosine(Rij, Rc, derivative=False):
    # Rij is the norm
    if derivative is False:
        result = 0.5 * (np.cos(np.pi * Rij / Rc) + 1.)
    else:
        result = -0.5 * np.pi / Rc * np.sin(np.pi * Rij / Rc)
    return result


class MyDescriptorSO3:
    '''
    A class to generate the SO3 power spectrum components
    based off of the Gaussian atomic neighbor density function
    defined in "On Representing Atomic Environments".

    args:
        nmax: int, degree of radial expansion
        lmax: int, degree of spherical harmonic expansion
        rcut: float, cutoff radius for neighbor calculation
        alpha: float, gaussian width parameter
        derivative: bool, whether to calculate the gradient of not
        weight_on: bool, if True, the neighbors with different type will be counted as negative
        primitive: bool, use the asePrimitiveNeighborList
    '''

    def __init__(self,
                 nmax=3, lmax=3, rcut=3.5, alpha=2.0,
                 derivative=True, stress=False,
                 cutoff_function='cosine', weight_on=False, primitive=False):
        # populate attributes
        self.nmax = nmax
        self.lmax = lmax
        self.rcut = rcut
        self.alpha = alpha
        self.derivative = derivative
        self.stress = stress
        self._type = "SO3"
        self.cutoff_function = cutoff_function
        self.weight_on = weight_on
        self.primitive = primitive
        #return

    def __str__(self):
        s = "SO3 descriptor with Cutoff: {:6.3f}".format(self.rcut)
        s += " lmax: {:d}, nmax: {:d}, alpha: {:.3f}\n".format(self.lmax, self.nmax, self.alpha)
        return s

    @property
    def rcut(self):
        return self._rcut


    @rcut.setter
    def rcut(self, rcut):
        if isinstance(rcut, float) is True or isinstance(rcut, int) is True:
            if rcut <= 0:
                raise ValueError('rcut must be greater than zero')
            self._rcut = rcut
        else:
            raise ValueError('rcut must be a float')


    @property
    def nmax(self):
        return self._nmax

    @nmax.setter
    def nmax(self, nmax):
        if isinstance(nmax, int) is True:
            if nmax < 1:
                raise ValueError('nmax must be greater than or equal to 1')
            if nmax > 11:
                raise ValueError('nmax > 11 yields complex eigenvalues which will mess up the calculation')
            self._nmax = nmax
        else:
            raise ValueError('nmax must be an integer')

    @property
    def stress(self):
        return self._stress

    @stress.setter
    def stress(self, stress):
        if isinstance(stress, bool) is True:
            self._stress = stress
        else:
            raise ValueError('stress must be a boolean value')

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, alpha):
        if isinstance(alpha, float) is True or isinstance(alpha, int) is True:
            if alpha <= 0:
                raise ValueError('alpha must be greater than zero')
            self._alpha = alpha
        else:
            raise ValueError('alpha must be a float')

    @property
    def lmax(self):
        return self._lmax

    @lmax.setter
    def lmax(self, lmax):
        if isinstance(lmax, int) is True:
            if lmax < 0:
                raise ValueError('lmax must be greater than or equal to zero')
            elif lmax > 32:
                raise NotImplementedError('''Currently we only support Wigner-D matrices and spherical harmonics
                for arguments up to l=32.  If you need higher functionality, raise an issue
                in our Github and we will expand the set of supported functions''')
            self._lmax = lmax
        else:
            raise ValueError('lmax must be an integer')

    @property
    def cutoff_function(self):
        return self._cutoff_function

    @cutoff_function.setter
    def cutoff_function(self, cutoff_function):
        if isinstance(cutoff_function, str) is True:
            # more conditions
            if cutoff_function == 'cosine':
                self._cutoff_function = Cosine
            else:
                raise NotImplementedError('The requested cutoff function has not been implemented')
        else:
            raise ValueError('You must specify the cutoff function as a string')

    def initialize_arrays(self):
        # number of atoms in periodic arrangement
        # for a crystal this will be the number of
        # atoms in the unit cell
        # for a cluster/molecule(s) this will be the total number
        # of atoms
        ncell = len(self._atoms) #self._atoms)
        # degree of spherical harmonic expansion
        lmax = self.lmax
        # degree of radial expansion
        nmax = self.nmax
        # number of unique power spectrum components
        # this is given by the triangular elements of
        # the radial expansion multiplied by the degree
        # of spherical harmonic expansion (including 0)
        ncoefs = nmax*(nmax+1)//2*(lmax+1)

        self._plist = np.zeros((ncell, ncoefs), dtype=np.float64)
        self._dplist = np.zeros((len(self.seq), ncoefs, 3), dtype=np.float64) #ffr: forces?
        self._pstress = np.zeros((len(self.seq), ncoefs, 3, 3), dtype=np.float64) #ffr: stresses?

        return


    def build_neighbor_list(self, atom_ids=None):
        '''
        Builds a neighborlist for the calculation of bispectrum components for
        a given ASE atoms object given in the calculate method.
        '''
        atoms = self._atoms
        if atom_ids is None:
            atom_ids = range(len(atoms))

        cutoffs = [self.rcut/2]*len(atoms)
        if self.primitive:
            print("Using PrimitiveNeighborList")
            nl = PrimitiveNeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
            nl.build(atoms.pbc, atoms.cell, atoms.get_scaled_positions())
        else:
            print("Using default NeighborList")
            nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
            nl.update(atoms)

        center_atoms = []
        neighbors = []
        neighbor_indices = []
        atomic_weights = []
        temp_indices = []

        for i in atom_ids:
            # get center atom position vector
            center_atom = atoms.positions[i]
            # get indices and cell offsets for each neighbor
            indices, offsets = nl.get_neighbors(i)
            temp_indices.append(indices)
            for j, offset in zip(indices, offsets):
                pos = atoms.positions[j] + np.dot(offset,atoms.get_cell()) - center_atom
                center_atoms.append(center_atom)
                neighbors.append(pos)
                if self.weight_on and atoms[j].number != atoms[i].number:
                    factor = -1
                else:
                    factor = 1
                atomic_weights.append(factor*atoms[j].number)
                neighbor_indices.append([i,j])

        neighbor_indices = np.array(neighbor_indices, dtype=np.int64)
        Seq = []
        for i in atom_ids:
            ineighs = neighbor_indices[:,0] == i
            unique_atoms = np.unique(neighbor_indices[ineighs])
            if i not in unique_atoms:
                at = list(unique_atoms)
                at.append(i)
                at.sort()
                unique_atoms = np.array(at)
            for j in unique_atoms:
                Seq.append([i,j])

        Seq = np.array(Seq, dtype=np.int64)
        self.center_atoms = np.array(center_atoms, dtype=np.float64)
        self.neighborlist = np.array(neighbors, dtype=np.float64)
        self.seq = Seq
        self.atomic_weights = np.array(atomic_weights, dtype=np.int64)
        self.neighbor_indices = neighbor_indices
        return


    def clear_memory(self):
        '''
        Clears all memory that isn't an essential attribute for the calculator
        '''
        attrs = list(vars(self).keys())
        for attr in attrs:
            if attr not in {'_nmax', '_lmax', '_rcut', '_alpha', '_derivative', '_stress', '_cutoff_function', 'weight_on', 'primitive'}:
                delattr(self, attr)
        return

    



def calculate_SO3_power_spectrum(SO3_calc, atoms, atom_ids=None):

    # Set a copy ?
    SO3_calc._atoms = atoms # ._atoms is not referenced anymore ?

    SO3_calc.build_neighbor_list(atom_ids)
    SO3_calc.initialize_arrays()

    # nmax*(nmax+1) // 2*(lmax+1)
    ncoefs = SO3_calc.nmax*(SO3_calc.nmax+1)//2*(SO3_calc.lmax+1)
    print("ncoefs = ", ncoefs)

    tril_indices = np.tril_indices(SO3_calc.nmax, k=0)

    ls = np.arange(SO3_calc.lmax+1)
    norm = np.sqrt(2*np.sqrt(2)*np.pi/np.sqrt(2*ls+1))

    if SO3_calc.derivative:
        # get expansion coefficients and derivatives
        cs, dcs = my_compute_dcs(
            SO3_calc.neighborlist,
            SO3_calc.nmax,
            SO3_calc.lmax,
            SO3_calc.rcut,
            SO3_calc.alpha,
            SO3_calc._cutoff_function # the this a function not str
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
        cs = my_compute_cs(
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
        x = {
            'x' : SO3_calc._plist,
            'dxdr' : None,
            'rdxdr' : None,
            'elements' : list(atoms.symbols)
        }

    #SO3_calc.clear_memory()
    return x


def calc_W(nmax):
    arr = np.zeros((nmax,nmax), np.float64)
    for alpha in range(1, nmax+1, 1):
        temp1 = (2*alpha+5)*(2*alpha+6)*(2*alpha+7)
        for beta in range(1, alpha+1, 1):
            temp2 = (2*beta+5)*(2*beta+6)*(2*beta+7)
            arr[alpha-1, beta-1] = np.sqrt(temp1*temp2)/(5+alpha+beta)/(6+alpha+beta)/(7+alpha+beta)
            arr[beta-1, alpha-1] = arr[alpha-1, beta-1]

    sinv = np.linalg.inv(arr)
    eigvals, V = np.linalg.eig(sinv)
    sqrtD = np.diag(np.sqrt(eigvals))
    arr[:,:] = np.dot(np.dot(V, sqrtD), np.linalg.inv(V))
    return arr


def GaussChebyshevQuadrature(nmax,lmax):
    NQuad = (nmax+lmax+1)*10
    quad_array = np.zeros(NQuad, dtype=np.float64)
    for i in range(1,NQuad+1,1):
        # roots of Chebyshev polynomial of degree N
        x = np.cos((2*i-1)*np.pi/2/NQuad)
        quad_array[i-1] = x
    return quad_array, np.pi/NQuad

# what's this?
def phi(r, alpha, rcut):
    return (rcut-r)**(alpha+2)/np.sqrt(2*rcut**(2*alpha+7)/(2*alpha+5)/(2*alpha+6)/(2*alpha+7))

# what's this?
def g(r, n, nmax, rcut, w):
    Sum = 0.0
    for alpha in range(1, nmax+1):
        Sum += w[n-1, alpha-1]*phi(r, alpha, rcut)
    return Sum


def my_compute_cs(pos, nmax, lmax, rcut, alpha, cutoff):

    print("Enter my_compute_cs")

    # compute the overlap matrix
    w = calc_W(nmax)

    # get the norm of the position vectors
    Ris = np.linalg.norm(pos, axis=1) # (Nneighbors)

    # initialize Gauss Chebyshev Quadrature
    GCQuadrature, weight = GaussChebyshevQuadrature(nmax,lmax) #(Nquad)
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
        Gs[n-1,:] = g(Quadrature, n, nmax, rcut, w)

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

    cutoff_array = cutoff(Ris, rcut)

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
            #ylms[:,l,midx] = sph_harm(m, l, phis, thetas) #ffr: old version
            ylms[:,l,midx] = sph_harm_y(l, m, thetas, phis)
    # In SciPy theta is the polar angle and phi is the azimuthal angle

    # multiply the spherical harmonics and the radial inner product
    Y_mul_innerprod = np.einsum('ijk,ilj->iljk', ylms, integral_array)

    # multiply the gaussians into the expression
    C = np.einsum('i,ijkl->ijkl', exparray, Y_mul_innerprod)
    return C


def my_compute_dcs(pos, nmax, lmax, rcut, alpha, cutoff):

    print("Enter my_compute_dcs")

    # compute the overlap matrix
    w = calc_W(nmax)
    print("w.shape = ", w.shape)

    # get the norm of the position vectors
    Ris = np.linalg.norm(pos, axis=1) # (Nneighbors)

    # get unit vectors
    upos = pos/Ris[:,np.newaxis]

    # initialize Gauss Chebyshev Quadrature
    GCQuadrature, weight = GaussChebyshevQuadrature(nmax,lmax) #(Nquad)
    weight *= rcut/2
    # transform from (-1,1) to (0, rcut)
    Quadrature = rcut/2*(GCQuadrature+1)

    # compute the arguments for the bessel functions
    BesselArgs = 2*alpha*np.outer(Ris,Quadrature)#(Nneighbors x Nquad)

    # initalize the arrays for the bessel function values
    # and the G function values
    Bessels = np.zeros((len(Ris), len(Quadrature), lmax+1), dtype=np.float64) #(Nneighbors x Nquad x lmax+1)
    Gs = np.zeros((nmax, len(Quadrature)), dtype=np.float64) # (nmax, nquad)
    dBessels = np.zeros((len(Ris), len(Quadrature), lmax+1), dtype=np.float64) #(Nneighbors x Nquad x lmax+1)

    # compute the g values
    for n in range(1,nmax+1,1):
        Gs[n-1,:] = g(Quadrature, n, nmax, rcut,w)*weight

    # compute the bessel values
    for l in range(lmax+1):
        Bessels[:,:,l] = spherical_in(l, BesselArgs)
        dBessels[:,:,l] = spherical_in(l, BesselArgs, derivative=True)

    #(Nneighbors x Nquad x lmax+1) unit vector here
    gradBessels = np.einsum('ijk,in->ijkn',dBessels,upos)
    gradBessels *= 2*alpha
    # multiply with r for the integral
    gradBessels = np.einsum('ijkn,j->ijkn',gradBessels,Quadrature)

    # mutliply the terms in the integral separate from the Bessels
    Quad_Squared = Quadrature**2
    Gs *= Quad_Squared * np.exp(-alpha*Quad_Squared) * np.sqrt(1-GCQuadrature**2)

    # perform the integration with the Bessels
    integral_array = np.einsum('ij,kjl->kil', Gs, Bessels) # (Nneighbors x nmax x lmax+1)

    grad_integral_array = np.einsum('ij,kjlm->kilm', Gs, gradBessels)# (Nneighbors x nmax x lmax+1, 3)

    # compute the gaussian for each atom
    exparray = 4*np.pi*np.exp(-alpha*Ris**2) # (Nneighbors)

    gradexparray = (-2*alpha*Ris*exparray)[:,np.newaxis]*upos

    cutoff_array = cutoff(Ris, rcut)

    grad_cutoff_array = np.einsum('i,in->in',cutoff(Ris, rcut, True), upos)

    # get the spherical coordinates of each atom
    thetas = np.arccos(pos[:,2]/Ris[:])
    phis = np.arctan2(pos[:,1], pos[:,0])

    # the size changes temporarily for the derivative
    # determine the size of the m axis
    Msize = 2*(lmax+1)+1
    msize = 2*lmax + 1
    # initialize an array for the spherical harmonics and gradients
    #(Nneighbors, l, m, *3*)
    ylms = np.zeros((len(Ris), lmax+1+1, Msize), dtype=np.complex128)
    gradylms = np.zeros((len(Ris), lmax+1, msize, 3), dtype=np.complex128)
    # compute the spherical harmonics
    for l in range(lmax+1+1):
        for m in range(-l,l+1,1):
            midx = Msize//2 + m
            ylms[:,l,midx] = sph_harm_y(l, m, thetas, phis)
    # thetas are polar angles

    for l in range(1, lmax+1):
        for m in range(-l, l+1, 1):
            midx = msize//2 + m
            Midx = Msize//2 + m
            # get gradient with recpect to spherical covariant components
            xcov0 = -np.sqrt(((l+1)**2-m**2)/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx]/Ris
            if abs(m) <= l-1:
                xcov0 += np.sqrt((l**2-m**2)/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx]/Ris
            #
            xcovpl1 = -np.sqrt((l+m+1)*(l+m+2)/2/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx+1]/Ris
            if abs(m+1) <= l-1:
                xcovpl1 -= np.sqrt((l-m-1)*(l-m)/2/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx+1]/Ris
            #
            xcovm1 = -np.sqrt((l-m+1)*(l-m+2)/2/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx-1]/Ris
            #
            if abs(m-1) <= l-1:
                xcovm1 -= np.sqrt((l+m-1)*(l+m)/2/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx-1]/Ris
            #
            #transform the gradient to cartesian
            gradylms[:,l,midx,0] = 1/np.sqrt(2)*(xcovm1-xcovpl1)
            gradylms[:,l,midx,1] = 1j/np.sqrt(2)*(xcovm1+xcovpl1)
            gradylms[:,l,midx,2] = xcov0

    # index ylms to get rid of extra terms for derivative
    ylms = ylms[:,0:lmax+1,1:1+2*lmax+1]
    # multiply the spherical harmonics and the radial inner product
    Y_mul_innerprod = np.einsum('ijk,ilj->iljk', ylms, integral_array)
    # multiply the gradient of the spherical harmonics with the radial inner get_radial_inner_product
    dY_mul_innerprod = np.einsum('ijkn,ilj->iljkn', gradylms, integral_array)
    # multiply the spherical harmonics with the gradient of the radial inner get_radial_inner_product
    Y_mul_dinnerprod = np.einsum('ijk,iljn->iljkn', ylms, grad_integral_array)
    # multiply the gaussians into the expression with 4pi
    C = np.einsum('i,ijkl->ijkl', exparray, Y_mul_innerprod)
    # multiply the gradient of the gaussian with the other terms
    gradexp_mul_y_inner = np.einsum('in,ijkl->ijkln', gradexparray, Y_mul_innerprod)
    # add gradient of inner product and spherical harmonic terms
    gradHarmonics_mul_gaussian = np.einsum('ijkln,i->ijkln', dY_mul_innerprod+Y_mul_dinnerprod, exparray)
    dC = gradexp_mul_y_inner + gradHarmonics_mul_gaussian
    dC *= cutoff_array[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
    dC += np.einsum('in,ijkl->ijkln', grad_cutoff_array, C)
    C *= cutoff_array[:,np.newaxis,np.newaxis,np.newaxis]
    return C, dC