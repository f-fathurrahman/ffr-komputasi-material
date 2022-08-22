# -*- coding: utf-8 -*-
import warnings

import numpy as np
from numpy.linalg import eigh

#from ase.optimize.optimize import Optimizer
from my_optimizer import MyOptimizer
from ase.utils import basestring

class MyBFGS(MyOptimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=0.04, master=None):
        """BFGS optimizer.

        Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: string
            Pickle file used to store hessian matrix. If set, file with
            such a name will be searched and hessian matrix stored will
            be used, if the file exists.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.04 Å).

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.
        """

        #print()
        #print("Initializing MyBFGS")
        #print("===================")

        if maxstep > 1.0:
            warnings.warn('You are using a much too large value for '
                          'the maximum step size: %.1f Å' % maxstep)
        self.maxstep = maxstep

        MyOptimizer.__init__(self, atoms, restart, logfile, trajectory, master)

    def todict(self):
        #print("MyBFGS: todict")
        d = MyOptimizer.todict(self)
        if hasattr(self, 'maxstep'):
            d.update(maxstep=self.maxstep)
        return d

    def initialize(self):
        #print("MyBFGS: initialize")
        self.H = None
        self.r0 = None
        self.f0 = None

    def read(self):
        #print("MyBFGS: read")
        self.H, self.r0, self.f0, self.maxstep = self.load()

    def step(self, f=None):
        print("MyBFGS: step")
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()
        #
        r = atoms.get_positions()
        f = f.reshape(-1)
        #
        self.update(r.flat, f, self.r0, self.f0)
        #
        omega, V = eigh(self.H)
        print("omega = ", omega)
        #
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        print("dr = ", dr)
        #
        steplengths = (dr**2).sum(1)**0.5
        print("MyBFGS: steplengths = ", steplengths)
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))
        #print("MyBFGS: end of step")

    def determine_step(self, dr, steplengths):
        """Determine step to take according to maxstep

        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        #print("MyBFGS: determine_step")
        maxsteplength = np.max(steplengths)
        if maxsteplength >= self.maxstep:
            dr *= self.maxstep / maxsteplength
        #print("MyBFGS: end of determine_step")
        return dr

    def update(self, r, f, r0, f0):
        #print("MyBFGS: update")
        if self.H is None:
            self.H = np.eye(3 * len(self.atoms)) * 70.0
            return
        dr = r - r0
        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        df = f - f0
        a = np.dot(dr, df)
        dg = np.dot(self.H, dr)
        b = np.dot(dr, dg)
        self.H -= np.outer(df, df) / a + np.outer(dg, dg) / b
        #print("MyBFGS: end of update")

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        #print("MyBFGS: replay_trajectory()")
        if isinstance(traj, basestring):
            from ase.io.trajectory import Trajectory
            traj = Trajectory(traj, 'r')
        self.H = None
        atoms = traj[0]
        r0 = atoms.get_positions().ravel()
        f0 = atoms.get_forces().ravel()
        for atoms in traj:
            r = atoms.get_positions().ravel()
            f = atoms.get_forces().ravel()
            self.update(r, f, r0, f0)
            r0 = r
            f0 = f

        self.r0 = r0
        self.f0 = f0


