from math import sqrt
import time

from ase.parallel import world, barrier
from my_dynamics import MyDynamics

class MyOptimizer(MyDynamics):
    """Base-class for all structure optimization classes."""

    def __init__(
        self,
        atoms,
        restart,
        logfile,
        trajectory,
        master=None,
        append_trajectory=False,
        force_consistent=False,
    ):
        """Structure optimizer object.

        Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: str
            Filename for restart file.  Default value is *None*.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        trajectory: Trajectory object or str
            Attach trajectory object.  If *trajectory* is a string a
            Trajectory will be constructed.  Use *None* for no
            trajectory.

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.

        append_trajectory: boolean
            Appended to the trajectory file instead of overwriting it.

        force_consistent: boolean or None
            Use force-consistent energy calls (as opposed to the energy
            extrapolated to 0 K).  If force_consistent=None, uses
            force-consistent energies if available in the calculator, but
            falls back to force_consistent=False if not.
        """
        MyDynamics.__init__(
            self,
            atoms,
            logfile,
            trajectory,
            append_trajectory=append_trajectory,
            master=master,
        )

        self.force_consistent = force_consistent
        if self.force_consistent is None:
            self.set_force_consistent()

        self.restart = restart

        # initialize attribute
        self.fmax = None

        if restart is None or not isfile(restart):
            self.initialize()
        else:
            self.read()
            barrier()

    def todict(self):
        description = {
            "type": "optimization",
            "optimizer": self.__class__.__name__,
        }
        return description

    def initialize(self):
        pass

    def irun(self, fmax=0.05, steps=None):
        """ call Dynamics.irun and keep track of fmax"""
        self.fmax = fmax
        if steps:
            self.max_steps = steps
        return MyDynamics.irun(self)

    def run(self, fmax=0.05, steps=None):
        """ call Dynamics.run and keep track of fmax"""
        self.fmax = fmax
        if steps:
            self.max_steps = steps
        return MyDynamics.run(self)

    def converged(self, forces=None):
        #print("MyOptimizer: converged?")
        """Did the optimization converge?"""
        if forces is None:
            #print("MyOptimizer: Calculating forces (again?)")
            forces = self.atoms.get_forces()
        if hasattr(self.atoms, "get_curvature"):
            return (forces ** 2).sum(
                axis=1
            ).max() < self.fmax ** 2 and self.atoms.get_curvature() < 0.0
        #print("MyOptimizer: end of converged?")
        max_forces = (forces ** 2).sum(axis=1).max()
        print("Check converged: max_forces = %18.10f, fmax = %18.10f" % (max_forces, self.fmax**2))
        #
        return max_forces < self.fmax ** 2

    def log(self, forces=None):
        #print("MyOptimizer: log")
        if forces is None:
            #print("MyOptimizer: calculating force (again?)")
            forces = self.atoms.get_forces()
        fmax = sqrt((forces ** 2).sum(axis=1).max())
        e = self.atoms.get_potential_energy(
            force_consistent=self.force_consistent
        )
        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                args = (" " * len(name), "Step", "Time", "Energy", "fmax")
                msg = "%s  %4s %8s %15s %12s\n" % args
                self.logfile.write(msg)

                if self.force_consistent:
                    msg = "*Force-consistent energies used in optimization.\n"
                    self.logfile.write(msg)

            ast = {1: "*", 0: ""}[self.force_consistent]
            args = (name, self.nsteps, T[3], T[4], T[5], e, ast, fmax)
            msg = "%s:  %3d %02d:%02d:%02d %15.6f%1s %12.4f\n" % args
            self.logfile.write(msg)

            self.logfile.flush()
        #print("MyOptimizer: end of log")

    def dump(self, data):
        if world.rank == 0 and self.restart is not None:
            pickle.dump(data, open(self.restart, "wb"), protocol=2)

    def load(self):
        return pickle.load(open(self.restart, "rb"))

    def set_force_consistent(self):
        """Automatically sets force_consistent to True if force_consistent
        energies are supported by calculator; else False."""
        try:
            self.atoms.get_potential_energy(force_consistent=True)
        except PropertyNotImplementedError:
            self.force_consistent = False
        else:
            self.force_consistent = True
