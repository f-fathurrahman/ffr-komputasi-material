"""Defines the System Class, derived from base.system module.

The System is supposed to be a singleton that
is passed to a singleton NeighborKMC object.

See Also  
---------
Module: base.system
Module: user_sites

"""

import numpy as np
from base.system import SystemBase
from ase import neighborlist


class System(SystemBase):
    """Class defines a collection of sites and connected atoms.
            
    Calls the base class system.py constructor, 
    sets the global neighborlist from the individual site's
    neighborlist.

    Attributes
    -----------
    atoms: ase.Atoms
        Can (optionally) be passed to connect an ASE.Atoms
        object to the system. This can be useful for visualization
        of the simulation, for example by setting the ase.Atoms tag
        according to coverages.
    sites: list(Site)
        The sites that constitute the system.

    See Also
    ---------
    Module: base.system

    """

    def __init__(self, atoms=None, sites=[]):
        SystemBase.__init__(self, atoms=atoms, sites=sites)

    def set_neighbors(self,  pbc=True):
        """Sets neighborlists of self.sites by distances.

        Loops through the sites and using self.atoms, the
        method adds neighbors to the sites that are within a
        neighbor-distance (Ncutoff).

        Parameters
        -----------
        pbc: bool
            If the neighborlist should be computed with periodic boundary
            conditions. To make a direction aperiodic, introduce a vacuum larger
            than Ncutoff in this direction in self.atoms.

        Raises
        ---------
        Warning:
            If self.atoms is not set, because then distances cannot
            be used to determine neighbors.

        """
        if self.atoms is None:
            raise Warning("Tried to set neighbor-distances in user_system.set_neighbors() with self.atom = None")
        
        cutOff = neighborlist.natural_cutoffs(self.atoms)
        nl = neighborlist.NeighborList(cutOff, self_interaction = False, bothways = True)
        nl.update(self.atoms)

        # Set the neighbor list for each site using distances.
        # ------------------------------------------
        for i, s in enumerate(self.sites):
            neigh_list = list(nl.get_neighbors(s.ind)[0])
            for j in neigh_list:
                s.neighbors.append(j)

    def cover_system(self, species, coverage):
        """Covers the system with a certain species.
            
        Randomly covers the system with a species, at a
        certain fractional coverage.
    
        Parameters
        ----------
        species: int
            The species as defined by the user (e.g. empty=0,CO=1).
        coverage: float
            The fractional coverage to load lattice with.

        """
        n_covered = int(np.round(coverage * len(self.system.sites)))
        chosen_sites = np.random.choice(len(self.system.sites), n_covered)
        for c in chosen_sites:
            self.system.sites[c].covered = species
