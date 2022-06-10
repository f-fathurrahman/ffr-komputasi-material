
TightBindingModel
-----------------

__init__(self,dim_k,dim_r,lat=None,orb=None,per=None,nspin=1)

This is the main class of the PythTB package which contains all
information for the tight-binding model.

:param dim_k: Dimensionality of reciprocal space, i.e., specifies how
  many directions are considered to be periodic.

:param dim_r: Dimensionality of real space, i.e., specifies how many
  real space lattice vectors there are and how many coordinates are
  needed to specify the orbital coordinates.
.. note:: Parameter *dim_r* can be larger than *dim_k*! For example,
  a polymer is a three-dimensional molecule (one needs three
  coordinates to specify orbital positions), but it is periodic
  along only one direction. For a polymer, therefore, we should
  have *dim_k* equal to 1 and *dim_r* equal to 3. See similar example
  here: :ref:`trestle-example`.

:param lat: Array containing lattice vectors in Cartesian
  coordinates (in arbitrary units). In example the below, the first
  lattice vector has coordinates [1.0,0.5] while the second
  one has coordinates [0.0,2.0].  By default, lattice vectors
  are an identity matrix.

:param orb: Array containing reduced coordinates of all
  tight-binding orbitals. In the example below, the first
  orbital is defined with reduced coordinates [0.2,0.3]. Its
  Cartesian coordinates are therefore 0.2 times the first
  lattice vector plus 0.3 times the second lattice vector.
  If *orb* is an integer code will assume that there are these many
  orbitals all at the origin of the unit cell.  By default
  the code will assume a single orbital at the origin.

:param per: This is an optional parameter giving a list of lattice
  vectors which are considered to be periodic. In the example below,
  only the vector [0.0,2.0] is considered to be periodic (since
  per=[1]). By default, all lattice vectors are assumed to be
  periodic. If dim_k is smaller than dim_r, then by default the first
  dim_k vectors are considered to be periodic.

:param nspin: Number of explicit spin components assumed for each
  orbital in *orb*. Allowed values of *nspin* are *1* and *2*. If
  *nspin* is 1 then the model is spinless, if *nspin* is 2 then it
  is explicitly a spinfull model and each orbital is assumed to
  have two spin components. Default value of this parameter is
  *1*.  Of course one can make spinfull calculation even with
  *nspin* set to 1, but then the user must keep track of which
  orbital corresponds to which spin component.

Example usage::

  # Creates model that is two-dimensional in real space but only
  # one-dimensional in reciprocal space. Second lattice vector is
  # chosen to be periodic (since per=[1]). Three orbital
  # coordinates are specified.       
  tb = TightBindingModel(1, 2,
              lat=[[1.0, 0.5], [0.0, 2.0]],
              orb=[[0.2, 0.3], [0.1, 0.1], [0.2, 0.2]],
              per=[1])





TightBindingModel.set_onsite
----------------------------

set_onsite(self, onsite_en, ind_i=None, mode="set"):

Defines on-site energies for tight-binding orbitals. One can
either set energy for one tight-binding orbital, or all at
once.

.. warning:: In previous version of PythTB this function was
  called *set_sites*. For backwards compatibility one can still
  use that name but that feature will be removed in future
  releases.

:param onsite_en: Either a list of on-site energies (in
  arbitrary units) for each orbital, or a single on-site
  energy (in this case *ind_i* parameter must be given). In
  the case when *nspin* is *1* (spinless) then each on-site
  energy is a single number.  If *nspin* is *2* then on-site
  energy can be given either as a single number, or as an
  array of four numbers, or 2x2 matrix. If a single number is
  given, it is interpreted as on-site energy for both up and
  down spin component. If an array of four numbers is given,
  these are the coefficients of I, sigma_x, sigma_y, and
  sigma_z (that is, the 2x2 identity and the three Pauli spin
  matrices) respectively. Finally, full 2x2 matrix can be
  given as well. If this function is never called, on-site
  energy is assumed to be zero.

:param ind_i: Index of tight-binding orbital whose on-site
  energy you wish to change. This parameter should be
  specified only when *onsite_en* is a single number (not a
  list).
          
:param mode: Similar to parameter *mode* in function set_hop*. 
  Speficies way in which parameter *onsite_en* is
  used. It can either set value of on-site energy from scratch,
  reset it, or add to it.

  * "set" -- Default value. On-site energy is set to value of
    *onsite_en* parameter. One can use "set" on each
    tight-binding orbital only once.

  * "reset" -- Specifies on-site energy to given value. This
    function can be called multiple times for the same
    orbital(s).

  * "add" -- Adds to the previous value of on-site
    energy. This function can be called multiple times for the
    same orbital(s).

Example usage::

  # Defines on-site energy of first orbital to be 0.0,
  # second 1.0, and third 2.0
  tb.set_onsite([0.0, 1.0, 2.0])
  # Increases value of on-site energy for second orbital
  tb.set_onsite(100.0, 1, mode="add")
  # Changes on-site energy of second orbital to zero
  tb.set_onsite(0.0, 1, mode="reset")
  # Sets all three on-site energies at once
  tb.set_onsite([2.0, 3.0, 4.0], mode="reset")



TightBinding.set_hop
--------------------

set_hop(self, hop_amp, ind_i, ind_j, ind_R=None, mode="set", allow_conjugate_pair=False):
        
Defines hopping parameters between tight-binding orbitals. In
the notation used in section 3.1 equation 3.6 of
:download:`notes on tight-binding formalism
<misc/pythtb-formalism.pdf>` this function specifies the
following object

.. math::

  H_{ij}({\bf R})= \langle \phi_{{\bf 0} i}  \vert H  \vert \phi_{{\bf R},j} \rangle

Where :math:`\langle \phi_{{\bf 0} i} \vert` is i-th
tight-binding orbital in the home unit cell and
:math:`\vert \phi_{{\bf R},j} \rangle` is j-th tight-binding orbital in
unit cell shifted by lattice vector :math:`{\bf R}`. :math:`H`
is the Hamiltonian.

(Strictly speaking, this term specifies hopping amplitude
for hopping from site *j+R* to site *i*, not vice-versa.)

Hopping in the opposite direction is automatically included by
the code since

.. math::

  H_{ji}(-{\bf R})= \left[ H_{ij}({\bf R}) \right]^{*}

.. warning::

  There is no need to specify hoppings in both :math:`i
  \rightarrow j+R` direction and opposite :math:`j
  \rightarrow i-R` direction since that is done
  automatically. If you want to specifiy hoppings in both
  directions, see description of parameter
  *allow_conjugate_pair*.

.. warning:: In previous version of PythTB this function was
  called *add_hop*. For backwards compatibility one can still
  use that name but that feature will be removed in future
  releases.

:param hop_amp: Hopping amplitude; can be real or complex
  number, equals :math:`H_{ij}({\bf R})`. If *nspin* is *2*
  then hopping amplitude can be given either as a single
  number, or as an array of four numbers, or as 2x2 matrix. If
  a single number is given, it is interpreted as hopping
  amplitude for both up and down spin component.  If an array
  of four numbers is given, these are the coefficients of I,
  sigma_x, sigma_y, and sigma_z (that is, the 2x2 identity and
  the three Pauli spin matrices) respectively. Finally, full
  2x2 matrix can be given as well.

:param ind_i: Index of bra orbital from the bracket :math:`\langle
  \phi_{{\bf 0} i} \vert H \vert \phi_{{\bf R},j} \rangle`. This
  orbital is assumed to be in the home unit cell.

:param ind_j: Index of ket orbital from the bracket :math:`\langle
  \phi_{{\bf 0} i} \vert H \vert \phi_{{\bf R},j} \rangle`. This
  orbital does not have to be in the home unit cell; its unit cell
  position is determined by parameter *ind_R*.

:param ind_R: Specifies, in reduced coordinates, the shift of
  the ket orbital. The number of coordinates must equal the
  dimensionality in real space (*dim_r* parameter) for consistency,
  but only the periodic directions of ind_R will be considered. If
  reciprocal space is zero-dimensional (as in a molecule),
  this parameter does not need to be specified.

:param mode: Similar to parameter *mode* in function *set_onsite*. 
  Speficies way in which parameter *hop_amp* is
  used. It can either set value of hopping term from scratch,
  reset it, or add to it.

  * "set" -- Default value. Hopping term is set to value of
    *hop_amp* parameter. One can use "set" for each triplet of
    *ind_i*, *ind_j*, *ind_R* only once.

  * "reset" -- Specifies on-site energy to given value. This
    function can be called multiple times for the same triplet
    *ind_i*, *ind_j*, *ind_R*.

  * "add" -- Adds to the previous value of hopping term This
    function can be called multiple times for the same triplet
    *ind_i*, *ind_j*, *ind_R*.

  If *set_hop* was ever called with *allow_conjugate_pair* set
  to True, then it is possible that user has specified both
  :math:`i \rightarrow j+R` and conjugate pair :math:`j
  \rightarrow i-R`.  In this case, "set", "reset", and "add"
  parameters will treat triplet *ind_i*, *ind_j*, *ind_R* and
  conjugate triplet *ind_j*, *ind_i*, *-ind_R* as distinct.

:param allow_conjugate_pair: Default value is *False*. If set
  to *True* code will allow user to specify hopping
  :math:`i \rightarrow j+R` even if conjugate-pair hopping
  :math:`j \rightarrow i-R` has been
  specified. If both terms are specified, code will
  still count each term two times.
          

Example usage::

  # Specifies complex hopping amplitude between first orbital in home
  # unit cell and third orbital in neigbouring unit cell.
  tb.set_hop(0.3+0.4j, 0, 2, [0, 1])
  # change value of this hopping
  tb.set_hop(0.1+0.2j, 0, 2, [0, 1], mode="reset")
  # add to previous value (after this function call below,
  # hopping term amplitude is 100.1+0.2j)
  tb.set_hop(100.0, 0, 2, [0, 1], mode="add")



TightBindingModel.k_path
------------------------

def k_path(self, kpts, nk, report=True)

Interpolates a path in reciprocal space between specified
k-points.  In 2D or 3D the k-path can consist of several
straight segments connecting high-symmetry points ("nodes"),
and the results can be used to plot the bands along this path.
        
The interpolated path that is returned contains as
equidistant k-points as possible.
    
:param kpts: Array of k-vectors in reciprocal space between
  which interpolated path should be constructed. These
  k-vectors must be given in reduced coordinates.  As a
  special case, in 1D k-space kpts may be a string:
    
  * *"full"*  -- Implies  *[ 0.0, 0.5, 1.0]*  (full BZ)
  * *"fullc"* -- Implies  *[-0.5, 0.0, 0.5]*  (full BZ, centered)
  * *"half"*  -- Implies  *[ 0.0, 0.5]*  (half BZ)
    
:param nk: Total number of k-points to be used in making the plot.
        
:param report: Optional parameter specifying whether printout
  is desired (default is True).

:returns:

  * **k_vec** -- Array of (nearly) equidistant interpolated
    k-points. The distance between the points is calculated in
    the Cartesian frame, however coordinates themselves are
    given in dimensionless reduced coordinates!  This is done
    so that this array can be directly passed to function
    :func:`pythtb.TightBindingModel.solve_all`.

  * **k_dist** -- Array giving accumulated k-distance to each
    k-point in the path.  Unlike array *k_vec* this one has
    dimensions! (Units are defined here so that for an
    one-dimensional crystal with lattice constant equal to for
    example *10* the length of the Brillouin zone would equal
    *1/10=0.1*.  In other words factors of :math:`2\pi` are
    absorbed into *k*.) This array can be used to plot path in
    the k-space so that the distances between the k-points in
    the plot are exact.

  * **k_node** -- Array giving accumulated k-distance to each
    node on the path in Cartesian coordinates.  This array is
    typically used to plot nodes (typically special points) on
    the path in k-space.
    
Example usage::

  # Construct a path connecting four nodal points in k-space
  # Path will contain 401 k-points, roughly equally spaced
  path = [[0.0, 0.0], [0.0, 0.5], [0.5, 0.5], [0.0, 0.0]]
  (k_vec,k_dist,k_node) = my_model.k_path(path,401)
  # solve for eigenvalues on that path
  evals = tb.solve_all(k_vec)
  # then use evals, k_dist, and k_node to plot bandstructure
  # (see examples)



TightBindingMode.solve_all
--------------------------

solve_all(self,k_list=None,eig_vectors=False):

Solves for eigenvalues and (optionally) eigenvectors of the
tight-binding model on a given one-dimensional list of k-vectors.

.. note::

Eigenvectors (wavefunctions) returned by this
function and used throughout the code are exclusively given
in convention 1 as described in section 3.1 of
:download:`notes on tight-binding formalism
<misc/pythtb-formalism.pdf>`.  In other words, they
are in correspondence with cell-periodic functions
:math:`u_{n {\bf k}} ({\bf r})` not
:math:`\Psi_{n {\bf k}} ({\bf r})`.

.. note::

In some cases class :class:`pythtb.wf_array` provides a more
elegant way to deal with eigensolutions on a regular mesh of
k-vectors.

:param k_list: One-dimensional array of k-vectors. Each k-vector
  is given in reduced coordinates of the reciprocal space unit
  cell. For example, for real space unit cell vectors [1.0,0.0]
  and [0.0,2.0] and associated reciprocal space unit vectors
  [2.0*pi,0.0] and [0.0,pi], k-vector with reduced coordinates
  [0.25,0.25] corresponds to k-vector [0.5*pi,0.25*pi].
  Dimensionality of each vector must equal to the number of
  periodic directions (i.e. dimensionality of reciprocal space,
  *dim_k*).
  This parameter shouldn't be specified for system with
  zero-dimensional k-space (*dim_k* =0).

:param eig_vectors: Optional boolean parameter, specifying whether
  eigenvectors should be returned. If *eig_vectors* is True, then
  both eigenvalues and eigenvectors are returned, otherwise only
  eigenvalues are returned.

:returns:
  * **eval** -- Two dimensional array of eigenvalues for
    all bands for all kpoints. Format is eval[band,kpoint] where
    first index (band) corresponds to the electron band in
    question and second index (kpoint) corresponds to the k-point
    as listed in the input parameter *k_list*. Eigenvalues are
    sorted from smallest to largest at each k-point seperately.

    In the case when reciprocal space is zero-dimensional (as in a
    molecule) kpoint index is dropped and *eval* is of the format
    eval[band].

  * **evec** -- Three dimensional array of eigenvectors for
    all bands and all kpoints. If *nspin* equals 1 the format
    of *evec* is evec[band,kpoint,orbital] where "band" is the
    electron band in question, "kpoint" is index of k-vector
    as given in input parameter *k_list*. Finally, "orbital"
    refers to the tight-binding orbital basis function.
    Ordering of bands is the same as in *eval*.  
            
    Eigenvectors evec[n,k,j] correspond to :math:`C^{n {\bf
    k}}_{j}` from section 3.1 equation 3.5 and 3.7 of the
    :download:`notes on tight-binding formalism
    <misc/pythtb-formalism.pdf>`.

    In the case when reciprocal space is zero-dimensional (as in a
    molecule) kpoint index is dropped and *evec* is of the format
    evec[band,orbital].

    In the spinfull calculation (*nspin* equals 2) evec has
    additional component evec[...,spin] corresponding to the
    spin component of the wavefunction.

Example usage::

    # Returns eigenvalues for three k-vectors
    eval = tb.solve_all([[0.0, 0.0], [0.0, 0.2], [0.0, 0.5]])
    # Returns eigenvalues and eigenvectors for two k-vectors
    (eval, evec) = tb.solve_all([[0.0, 0.0], [0.0, 0.2]], eig_vectors=True)




TightBindingModel.visualize
---------------------------

visualize(self,
  dir_first,
  dir_second=None,
  eig_dr=None, draw_hoppings=True, ph_color="black"):

Rudimentary function for visualizing tight-binding model geometry,
hopping between tight-binding orbitals, and electron eigenstates.

If eigenvector is not drawn, then orbitals in home cell are drawn
as red circles, and those in neighboring cells are drawn with
different shade of red. Hopping term directions are drawn with
green lines connecting two orbitals. Origin of unit cell is
indicated with blue dot, while real space unit vectors are drawn
with blue lines.

If eigenvector is drawn, then electron eigenstate on each orbital
is drawn with a circle whose size is proportional to wavefunction
amplitude while its color depends on the phase. There are various
coloring schemes for the phase factor; see more details under
*ph_color* parameter. If eigenvector is drawn and coloring scheme
is "red-blue" or "wheel", all other elements of the picture are
drawn in gray or black.

:param dir_first: First index of Cartesian coordinates used for
  plotting.

:param dir_second: Second index of Cartesian coordinates used for
  plotting. For example if dir_first=0 and dir_second=2, and
  Cartesian coordinates of some orbital is [2.0,4.0,6.0] then it
  will be drawn at coordinate [2.0,6.0]. If dimensionality of real
  space (*dim_r*) is zero or one then dir_second should not be
  specified.

:param eig_dr: Optional parameter specifying eigenstate to
  plot. If specified, this should be one-dimensional array of
  complex numbers specifying wavefunction at each orbital in
  the tight-binding basis. If not specified, eigenstate is not
  drawn.

:param draw_hoppings: Optional parameter specifying whether to
  draw all allowed hopping terms in the tight-binding
  model. Default value is True.

:param ph_color: Optional parameter determining the way
  eigenvector phase factors are translated into color. Default
  value is "black". Convention of the wavefunction phase is as
  in convention 1 in section 3.1 of :download:`notes on
  tight-binding formalism  <misc/pythtb-formalism.pdf>`.  In
  other words, these wavefunction phases are in correspondence
  with cell-periodic functions :math:`u_{n {\bf k}} ({\bf r})`
  not :math:`\Psi_{n {\bf k}} ({\bf r})`.

  * "black" -- phase of eigenvectors are ignored and wavefunction
    is always colored in black.

  * "red-blue" -- zero phase is drawn red, while phases or pi or
    -pi are drawn blue. Phases in between are interpolated between
    red and blue. Some phase information is lost in this coloring
    becase phase of +phi and -phi have same color.

  * "wheel" -- each phase is given unique color. In steps of pi/3
    starting from 0, colors are assigned (in increasing hue) as:
    red, yellow, green, cyan, blue, magenta, red.

:returns:
  * **fig** -- Figure object from matplotlib.pyplot module
    that can be used to save the figure in PDF, EPS or similar
    format, for example using fig.savefig("name.pdf") command.
  * **ax** -- Axes object from matplotlib.pyplot module that can be
    used to tweak the plot, for example by adding a plot title
    ax.set_title("Title goes here").


Example usage::

  # Draws x-y projection of tight-binding model
  # tweaks figure and saves it as a PDF.
  (fig, ax) = tb.visualize(0, 1)
  ax.set_title("Title goes here")
  fig.savefig("model.pdf")

See also these examples: :ref:`edge-example`, :ref:`visualize-example`.



TightBinding.position_matrix()
------------------------------

def position_matrix(self, evec, dir)

Returns matrix elements of the position operator along
direction *dir* for eigenvectors *evec* at a single k-point.
Position operator is defined in reduced coordinates.

The returned object :math:`X` is

.. math::

  X_{m n {\bf k}}^{\alpha} = \langle u_{m {\bf k}} \vert
  r^{\alpha} \vert u_{n {\bf k}} \rangle

Here :math:`r^{\alpha}` is the position operator along direction
:math:`\alpha` that is selected by *dir*.

:param evec: Eigenvectors for which we are computing matrix
  elements of the position operator.  The shape of this array
  is evec[band,orbital] if *nspin* equals 1 and
  evec[band,orbital,spin] if *nspin* equals 2.

:param dir: Direction along which we are computing the center.
  This integer must not be one of the periodic directions
  since position operator matrix element in that case is not
  well defined.

:returns:
  * **pos_mat** -- Position operator matrix :math:`X_{m n}` as defined 
    above. This is a square matrix with size determined by number of bands
    given in *evec* input array.  First index of *pos_mat* corresponds to
    bra vector (*m*) and second index to ket (*n*).

Example usage::

  # diagonalizes Hamiltonian at some k-points
  (evals, evecs) = my_model.solve_all(k_vec,eig_vectors=True)
  # computes position operator matrix elements for 3-rd kpoint 
  # and bottom five bands along first coordinate
  pos_mat = my_model.position_matrix(evecs[:5,2], 0)

See also this example: :ref:`haldane_hwf-example`,


TightBinding.position_expectation
---------------------------------

position_expectation(self,evec,dir)

Returns diagonal matrix elements of the position operator.
These elements :math:`X_{n n}` can be interpreted as an
average position of n-th Bloch state *evec[n]* along
direction *dir*.  Generally speaking these centers are *not*
hybrid Wannier function centers (which are instead
returned by :func:`pythtb.TightBindingModel.position_hwf`).
        
See function :func:`pythtb.TightBindingModel.position_matrix` for
definition of matrix :math:`X`.

:param evec: Eigenvectors for which we are computing matrix
  elements of the position operator.  The shape of this array
  is evec[band,orbital] if *nspin* equals 1 and
  evec[band,orbital,spin] if *nspin* equals 2.

:param dir: Direction along which we are computing matrix
  elements.  This integer must not be one of the periodic
  directions since position operator matrix element in that
  case is not well defined.

:returns:
  * **pos_exp** -- Diagonal elements of the position operator matrix :math:`X`.
    Length of this vector is determined by number of bands given in *evec* input 
    array.

Example usage::

  # diagonalizes Hamiltonian at some k-points
  (evals, evecs) = my_model.solve_all(k_vec,eig_vectors=True)
  # computes average position for 3-rd kpoint 
  # and bottom five bands along first coordinate
  pos_exp = my_model.position_expectation(evecs[:5,2], 0)

See also this example: :ref:`haldane_hwf-example`.



TightBinding.position_hwf
-------------------------

position_hwf(self,evec,dir,hwf_evec=False,basis="orbital")

Returns eigenvalues and optionally eigenvectors of the
position operator matrix :math:`X` in either Bloch or orbital
basis.  These eigenvectors can be interpreted as linear
combinations of Bloch states *evec* that have minimal extent (or
spread :math:`\Omega` in the sense of maximally localized
Wannier functions) along direction *dir*. The eigenvalues are
average positions of these localized states. 

Note that these eigenvectors are not maximally localized
Wannier functions in the usual sense because they are
localized only along one direction.  They are also not the
average positions of the Bloch states *evec*, which are
instead computed by :func:`pythtb.TightBindingModel.position_expectation`.

See function :func:`pythtb.TightBindingModel.position_matrix` for
the definition of the matrix :math:`X`.

See also Fig. 3 in Phys. Rev. Lett. 102, 107603 (2009) for a
discussion of the hybrid Wannier function centers in the
context of a Chern insulator.

:param evec: Eigenvectors for which we are computing matrix
  elements of the position operator.  The shape of this array
  is evec[band,orbital] if *nspin* equals 1 and
  evec[band,orbital,spin] if *nspin* equals 2.

:param dir: Direction along which we are computing matrix
  elements.  This integer must not be one of the periodic
  directions since position operator matrix element in that
  case is not well defined.

:param hwf_evec: Optional boolean variable.  If set to *True* 
  this function will return not only eigenvalues but also 
  eigenvectors of :math:`X`. Default value is *False*.

:param basis: Optional parameter. If basis="bloch" then hybrid
  Wannier function *hwf_evec* is written in the Bloch basis.  I.e. 
  hwf[i,j] correspond to the weight of j-th Bloch state from *evec*
  in the i-th hybrid Wannier function.  If basis="orbital" and nspin=1 then
  hwf[i,orb] correspond to the weight of orb-th orbital in the i-th 
  hybrid Wannier function.  If basis="orbital" and nspin=2 then
  hwf[i,orb,spin] correspond to the weight of orb-th orbital, spin-th
  spin component in the i-th hybrid Wannier function.  Default value
  is "orbital".

:returns:
  * **hwfc** -- Eigenvalues of the position operator matrix :math:`X`
    (also called hybrid Wannier function centers).
    Length of this vector equals number of bands given in *evec* input 
    array.  Hybrid Wannier function centers are ordered in ascending order.
    Note that in general *n*-th hwfc does not correspond to *n*-th electronic
    state *evec*.

  * **hwf** -- Eigenvectors of the position operator matrix :math:`X`.
    (also called hybrid Wannier functions).  These are returned only if
    parameter *hwf_evec* is set to *True*.
    The shape of this array is [h,x] or [h,x,s] depending on value of *basis*
    and *nspin*.  If *basis* is "bloch" then x refers to indices of 
    Bloch states *evec*.  If *basis* is "orbital" then *x* (or *x* and *s*)
    correspond to orbital index (or orbital and spin index if *nspin* is 2).

Example usage::

  # diagonalizes Hamiltonian at some k-points
  (evals, evecs) = my_model.solve_all(k_vec,eig_vectors=True)
  # computes hybrid Wannier centers (and functions) for 3-rd kpoint 
  # and bottom five bands along first coordinate
  (hwfc, hwf) = my_model.position_hwf(evecs[:5,2], 0, hwf_evec=True, basis="orbital")

See also this example: :ref:`haldane_hwf-example`,
