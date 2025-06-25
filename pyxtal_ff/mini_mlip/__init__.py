import os
import torch

import numpy as np
from scipy.special import sph_harm, spherical_in

import ase.io
from ase.neighborlist import NeighborList, PrimitiveNeighborList

from copy import deepcopy

import gc
import sys
import time
import shelve
from torch import save, load

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use("ggplot")

eV2GPa = 160.21766

__version__ = "0.0.1"


"""
Module for handling atomic elements.
"""

class Element:
    """
    Class for storing and accessing atomic elements. 
    Args:
        input_value: The atomic number, symbol, or name of the element
    """
    def __init__(self, input_value):
        self.input = input_value

        # list with atomic number z, short name, full name, valence, 
                # valence electrons, covalent radius, vdW radius, metallic radius
        self.elements_list = [
            [1,   'H' ],  
            [2,   'He'],
            [3,   'Li'],
            [4,   'Be'],
            [5,   'B' ],
            [6,   'C' ],
            [7,   'N' ],
            [8,   'O' ],
            [9,   'F' ],
            [10,  'Ne'],
            [11,  'Na'],
            [12,  'Mg'],
            [13,  'Al'],
            [14,  'Si'],
            [15,  'P' ],
            [16,  'S' ],
            [17,  'Cl'],
            [18,  'Ar'],
            [19,  'K' ],
            [20,  'Ca'],
            [21,  'Sc'],
            [22,  'Ti'],
            [23,  'V' ],
            [24,  'Cr'],
            [25,  'Mn'],
            [26,  'Fe'],
            [27,  'Co'],
            [28,  'Ni'],
            [29,  'Cu'],
            [30,  'Zn'],
            [31,  'Ga'],
            [32,  'Ge'],
            [33,  'As'],
            [34,  'Se'],
            [35,  'Br'],
            [36,  'Kr'],
            [37,  'Rb'],
            [38,  'Sr'],
            [39,  'Y' ],
            [40,  'Zr'],
            [41,  'Nb'],
            [42,  'Mo'],
            [43,  'Tc'],
            [44,  'Ru'],
            [45,  'Rh'],
            [46,  'Pd'],
            [47,  'Ag'],
            [48,  'Cd'],
            [49,  'In'],
            [50,  'Sn'],
            [51,  'Sb'],
            [52,  'Te'],
            [53,  'I' ],
            [54,  'Xe'],
            [55,  'Cs'],
            [56,  'Ba'],
            [57,  'La'],
            [58,  'Ce'],
            [59,  'Pr'],
            [60,  'Nd'],
            [61,  'Pm'],
            [62,  'Sm'],
            [63,  'Eu'],
            [64,  'Gd'],
            [65,  'Tb'],
            [66,  'Dy'],
            [67,  'Ho'],
            [68,  'Er'],
            [69,  'Tm'],
            [70,  'Yb'],
            [71,  'Lu'],
            [72,  'Hf'],
            [73,  'Ta'],
            [74,  'W' ],
            [75,  'Re'],
            [76,  'Os'],
            [77,  'Ir'],
            [78,  'Pt'],
            [79,  'Au'],
            [80,  'Hg'],
            [81,  'Tl'],
            [82,  'Pb'],
            [83,  'Bi'],
            [84,  'Po'],
            [85,  'At'],
            [86,  'Rn'],
            [87,  'Fr'],
            [88,  'Ra'],
            [89,  'Ac'],
            [90,  'Th'],
            [91,  'Pa'],
            [92,  'U' ],
            [93,  'Np'],
            [94,  'Pu'],
            [95,  'Am'],
            [96,  'Cm'],
            [97,  'Bk'],
            [98,  'Cf'],
            [99,  'Es'],
            [100, 'Fm'],
            [101, 'Md'],
            [102, 'No'],
            [103, 'Lr'],
            [104, 'Rf'],
            [105, 'Db'],
        ]
        """A list of atomic numbers, symbols, names, and other information, up
        to atomic number 105"""
        
        Z, self.Z = [], []
        for element in self.elements_list:
            if element[1] in self.input:
                Z.append(element)


        for inp in self.input:
            for z in Z:
                if inp == z[1]:
                    self.Z.append(z[0])

    def get_Z(self,):
        return self.Z






class MiniMLIP():

    def __init__(self, descriptors=None, model=None):
        
        # Checking the keys in descriptors
        descriptors_keywords = [
            "type", "Rc", "weights", "N_train", "N_test", "cutoff",
            "force", "stress", "ncpu", "parameters", "base_potential"
        ]
        
        if descriptors is not None:
            for key in descriptors.keys():
                if key not in descriptors_keywords:
                    msg = f"Don't recognize {key} in descriptors. "+\
                          f"Here are the keywords: {descriptors_keywords}."
                    raise NotImplementedError(msg)

        # Set up default descriptors parameters
        self._descriptors = {
            "system": model["system"],
            "type": "Bispectrum",
            "Rc": 5.0,
            "weights": None,
            "N": None,
            "N_train": None,
            "N_test": None,
            "ncpu": 1,
            "force": True,
            "stress": True,
            "cutoff": "cosine",
            "base_potential": False,
        }
        
        # descriptors is a Dict
        # Set default value
        if descriptors is not None:
            self._descriptors.update(descriptors)
            if "type" in descriptors and descriptors["type"] in ["EAD", "ead"]:
                _parameters = {
                    "L": 3,
                    "eta": [0.1],
                    "Rs": [1.]
                }
            elif "type" in descriptors and descriptors["type"] in ["SO3", "SOAP"]:
                _parameters = {
                    "nmax": 3,
                    "lmax": 3,
                    "alpha": 2.0
                }
            else:
                _parameters = {
                    "lmax": 3,
                    "rfac": 0.99363,
                    "normalize_U": False
                }
            
            if "parameters" in descriptors:
                _parameters.update(descriptors["parameters"])
                self._descriptors["parameters"] = _parameters

        # Create new directory to dump all the results.
        # E.g. for default "Si-O-Bispectrum/"
        if "path" in model:
            self.path = model["path"]
        else:
            _system = model["system"]
            self.path = "-".join(_system) + "-"
            self.path += self._descriptors["type"] + "/"

        if not os.path.exists(self.path):
            os.mkdir(self.path)

        self.print_descriptors(self._descriptors)
        
        # Checking the keys in model.
        keywords = ["algorithm", "system", "hiddenlayers", "activation", 
                    "random_seed", "force_coefficient", "unit", "softmax_beta", 
                    "restart", "optimizer", "path", "order", "d_max", 
                    "epoch", "device", "alpha", "batch_size", "noise", "kernel",
                    "norm", "stress_coefficient", "stress_group", "memory"]
        for key in model.keys():
            if key not in keywords:
                msg = f"Don't recognize {key} in model. "+\
                      f"Here are the keywords: {keywords}."
                raise NotImplementedError(msg)
        
        # Create model
        pr_keywords = ["PolynomialRegression", "PR"]
        nn_keywords = ["NeuralNetwork", "NN"]
        if "algorithm" not in model:
            model["algorithm"] = "NN"
        # default is using NN

        if model["algorithm"] in pr_keywords:
            self.algorithm = "PR"
        elif model["algorithm"] in nn_keywords:
            self.algorithm = "NN"
        else:
            msg = f"{model['algorithm']} is not implemented."
            raise NotImplementedError(msg)
        
        self._model = model

    def todict(self):
        return {"descriptor": self._descriptors, "model": self._model}



    def run(self, mode="train", TrainData=None, TestData=None, mliap=None):

        if mode == "train":
            assert TrainData is not None, "TrainData can't be None for train mode."

            # Instantiate model
            print("Initializing model ...")
            self._MODEL(self._model)
            print("... Done initializing model")

            # Calculate descriptors.
            self._descriptors.update({"N": self._descriptors["N_train"]})
            if not os.path.exists(self.path+"Train_db.dat") and not os.path.exists(self.path+"Train_db.db"):
                trainDB = Database(name=self.path+"Train_db")
                trainDB.store(TrainData, self._descriptors, True, self.path+"ase.db")
            else:
                # The database is already exist, just load it
                # XXX In several case this will error if the database is not complete
                trainDB = Database(name=self.path+"Train_db")
                trainDB.store(TrainData, self._descriptors, False)
            #
            trainDB.close()

            if TestData is not None:
                EvaluateTest = True
                self._descriptors.update({"N": self._descriptors["N_test"]}) 
                if not os.path.exists(self.path+"Test_db.dat"):
                    testDB = Database(name=self.path+"Test_db")
                    testDB.store(TestData, self._descriptors, True, self.path+"ase.db")
                else:
                    testDB = Database(name=self.path+"Test_db")
                    testDB.store(TestData, self._descriptors, False)
                testDB.close()

            else:
                EvaluateTest = False
            
            print("=========================== Training =============================\n")

            self.model.train("Train_db", optimizer=self.optimizer)
            self.model.save_checkpoint(des_info=self._descriptors)
            
            print("==================================================================\n")
            
            print(f"==================== Evaluating Training Set ====================\n")

            train_stat = self.model.evaluate("Train_db", figname="Train.png")
            
            print("==================================================================\n")

            if EvaluateTest:
                print("================= Evaluating Testing Set =====================\n")

                test_stat =  self.model.evaluate("Test_db", figname="Test.png")
                
                print("==============================================================\n")
            else:
                test_stat = None

            return (train_stat, test_stat)
        
        elif mode == "predict":
            #self._model["algorithm"] = torch.load(mliap)["algorithm"]
            # New in PyTorch 2.6 ?
            self._model["algorithm"] = torch.load(mliap, weights_only=False)["algorithm"]
            self.algorithm = self._model["algorithm"]
            self._MODEL(self._model)
            self._descriptors = self.model.load_checkpoint(filename=mliap)

    
    def _MODEL(self, model):
        """ Machine learning model is created here. """
                    
        if self.algorithm == "NN":
            raise NotImplementedError
                
        elif self.algorithm == "PR":
            _model = {
                "system": None,
                "force_coefficient": 0.0001,
                "stress_coefficient": None,
                "stress_group": None,
                "order": 1,
                "path": self.path,
                "alpha": None,
                "norm": 2,
                "d_max": None,
            }
            _model.update(model)
            self.model = PolynomialRegression(
                elements=_model["system"],
                force_coefficient=_model["force_coefficient"],
                stress_coefficient=_model["stress_coefficient"],
                stress_group=_model["stress_group"],
                order=_model["order"],
                path=_model["path"],
                alpha=_model["alpha"],
                norm=_model["norm"],
                d_max=_model["d_max"]
            )
            self.optimizer = None

    def print_descriptors(self, _descriptors):
        """ Print the descriptors information. """

        print("Descriptor parameters:")
        keys = ["type", "Rc", "cutoff"]
        for key in keys:
            print("{:12s}: {:}".format(key, _descriptors[key]))

        if _descriptors["type"] in ["SO4", "Bispectrum"]:
            key_params = ["lmax", "normalize_U"]
        elif _descriptors["type"] in ["SO3", "SOAP"]:
            key_params = ["nmax", "lmax", "alpha"]
        elif _descriptors["type"] in ["SNAP", "snap"]:
            key_params = ["lmax", "rfac"]
        elif _descriptors["type"] == "EAD":
            key_params = ["L", "eta", "Rs"]
        else:
            key_params = []

        for key in key_params:
            print("{:12s}: {:}".format(key, _descriptors["parameters"][key]))
        print("\n")





#import numpy as np
#from ase.optimize import LBFGS
#from ase.optimize.fire import FIRE
#from ase.constraints import ExpCellFilter
#from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
#np.set_printoptions(formatter={'float': '{: 8.4f}'.format})

from ase import units
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixSymmetry
from ase.filters import ExpCellFilter
from ase.optimize.fire import FIRE

class MiniMLIPCalculator(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']
    nolabel = True

    def __init__(self, style='ase', **kwargs):
        self.style = style
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):

        Calculator.calculate(self, atoms, properties, system_changes)

        #chem_symbols = list(set(atoms.get_chemical_symbols()))
        #self.ff = PyXtal_FF(model={'system': chem_symbols}, logo=self.parameters.logo)
        #self.ff.run(mode='predict', mliap=self.parameters.mliap)

        # base potential
        if self.parameters.ff._descriptors['base_potential']:
            raise NotImplementedError("Not implemented")
            #self.base_potential = ZBL(self.parameters.ff._descriptors['base_potential']['inner'],
            #                          self.parameters.ff._descriptors['base_potential']['outer'],
            #                          atomic_energy=True)
            #base_results = self.base_potential.calculate(atoms)
            #base_energy = base_results['energy']
            #base_forces = base_results['force']
            #base_stress = base_results['stress'] # eV/A^3
            #base_energies = base_results['energies']
        else:
            base_energy = 0
            base_forces = np.zeros([len(atoms), 3])
            base_stress = np.zeros([6])
            base_energies = 0.

        desp = compute_descriptor(self.parameters.ff._descriptors, atoms)
        energies, forces, stress = self.parameters.ff.model.calculate_properties(desp, bforce=True, bstress=True)

        self.desp = desp
        self.results['energies'] = energies + base_energies
        self.results['energy'] = energies.sum() + base_energy
        self.results['free_energy'] = energies.sum() + base_energy
        self.results['forces'] = forces + base_forces

        # pyxtal_ff and lammps uses: xx, yy, zz, xy, xz, yz
        # ase uses: xx, yy, zz, yz, xz, xy
        # vasp uses: xx, yy, zz, xy, yz, zx
        # from eV/A^3 to GPa 
        self.results['stress_zbl'] = base_stress/ase.units.GPa
        self.results['energy_zbl'] = base_energy
        self.results['forces_zbl'] = base_forces
        self.results['stress_ml'] = stress 
        self.results['energy_ml'] = energies.sum()
        self.results['forces_ml'] = forces


        # ase counts the stress differently
        if self.style == 'ase':
            self.results['stress'] = -(stress * units.GPa + base_stress)[[0, 1, 2, 5, 4, 3]]
        else:
            self.results['stress'] = self.results['stress_zbl'] + self.results['stress_ml']

    def __str__(self):
        s = "\nASE calculator with pyxtal_ff force field\n"
        return s

    def __repr__(self):
        return str(self)

    def print_stresses(self):
        print("stress_ml (GPa, xx, yy, zz, xy, xz, yz):", self.results["stress_ml"])
        print("stress_zbl(GPa, xx, yy, zz, xy, xz, yz):", self.results['stress_zbl'])

    def print_energy(self):
        print("energy_ml (eV):", self.results["energy_ml"])
        print("energy_zbl(eV):", self.results['energy_zbl'])

    def print_forces(self):
        print("forces (eV/A)")
        for f1, f2 in zip(self.results["forces_ml"], self.results['forces_zbl']):
            print("{:8.3f} {:8.3f} {:8.3f} -> {:8.3f} {:8.3f} {:8.3f}".format(*f1, *f2))

    def print_all(self):
        self.print_energy()
        self.print_forces()
        self.print_stresses()


def mini_mlip_elastic_properties(C):
    Kv = C[:3,:3].mean()
    Gv = (C[0,0]+C[1,1]+C[2,2] - (C[0,1]+C[1,2]+C[2,0]) + 3*(C[3,3]+C[4,4]+C[5,5]))/15
    Ev = 1/((1/(3*Gv))+(1/(9*Kv)))
    vv  = 0.5*(1-((3*Gv)/(3*Kv+Gv))); 

    S = np.linalg.inv(C)
    Kr = 1/((S[0,0]+S[1,1]+S[2,2])+2*(S[0,1]+S[1,2]+S[2,0])) 
    Gr = 15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[1,2]+S[2,0])+3*(S[3,3]+S[4,4]+S[5,5])) 
    Er = 1/((1/(3*Gr))+(1/(9*Kr))) 
    vr = 0.5*(1-((3*Gr)/(3*Kr+Gr))) 

    Kh = (Kv+Kr)/2    
    Gh = (Gv+Gr)/2    
    Eh = (Ev+Er)/2    
    vh = (vv+vr)/2   
    return Kv, Gv, Ev, vv, Kr, Gr, Er, vr, Kh, Gh, Eh, vh

def mini_mlip_optimize(atoms, sym=True, box=False, P=0.0, method='FIRE', fmax=0.01, steps=1000, logfile='ase.log'):
    """
    Geometry relaxation

    Args:
        Atoms: ase atoms
        sym: whether or not fix symmetry
        box: whether or not relax box
        P: external pressure in GPa
        method: optimization method
        fmax: toleration force
        steps: maximum number of steps
        logfile: output of the log file
    """
    if sym:
        atoms.set_constraint(FixSymmetry(atoms))
    if box:
        ecf = ExpCellFilter(atoms, scalar_pressure=P*units.GPa)
        if method == 'FIRE':
            dyn = FIRE(ecf, logfile=logfile)
        else:
            dyn = LBFGS(ecf, logfile=logfile)
    else:
        if method == 'FIRE':
            dyn = FIRE(atoms, logfile=logfile)
        else:
            dyn = FIRE(atoms, logfile=logfile)

    dyn.run(fmax=fmax, steps=steps)
    atoms.set_constraint()
    return atoms









class PolynomialRegression():
    """ Atom-centered Polynomial Regression (PR) model. PR utilizes
    linear regression to predict the energy and forces based on 
    the atom-centered descriptors as the input values.
    """
    def __init__(self,
            elements, force_coefficient, stress_coefficient, 
            stress_group, order, path, alpha, norm, d_max=None):
        
        self.force_coefficient = force_coefficient
        self.stress_coefficient = stress_coefficient
        self.stress_group = stress_group
        self.elements = sorted(elements)
        self.order = order
        if order == 1:
            self.quadratic = False
        elif order == 2:
            self.quadratic = True
        else:
            raise ValueError("Order must be 1 or 2")
        self.path = path

        self.alpha = alpha
        self.norm = norm
        self.d_max = d_max
        self.unit = 'eV'

    def train(self, TrainData, optimizer):
        """ Fitting Linear Regression model. """
        db = shelve.open(self.path+TrainData)
        self.no_of_structures = len(list(db.keys()))

        # d_max is the total number of descriptors used.
        if self.d_max is None:
            self.d_max = db['0']['x'].shape[1]
        else:
            # d_max has to be less or equal than total descriptors.
            assert self.d_max <= len(db['0']['x'].shape[1]),\
                    "d_max is larger than total descriptors."

        if self.stress_coefficient and (self.stress_group is None):
            sg = []
            for i in range(self.no_of_structures):
                if db[str(i)]['group'] not in sg:
                    sg.append(db[str(i)]['group'])
            self.stress_group = sg

        db.close()
        
        print(f"Order              : {self.order}")
        if self.order == 1:
            print(f"No of parameters   : {self.d_max+1}")
        else:
            print(f"No of parameters   : {(self.d_max+1)*(self.d_max+2)//2}")
        print(f"No of structures   : {self.no_of_structures}")
        print(f"Force_coefficient  : {self.force_coefficient}")
        print(f"Stress_coefficient : {self.stress_coefficient}")
        print(f"alpha              : {self.alpha}")
        print(f"norm               : {self.norm}\n")

        t0 = time.time()
        y, w = self.parse_features(TrainData)

        X = self.parse_descriptors(TrainData,
            fc=self.force_coefficient, sc=self.stress_coefficient)
        
        self.coef_ = self.LinearRegression(X, y, w, self.alpha, self.norm)
        
        t1 = time.time()
        print("The training time: {:.2f} s".format(t1-t0))
        
    
    def evaluate(self, data, figname):
        """ Evaluating the train or test data set. """
        db = shelve.open(self.path+data)
        
        energy, force, stress = [], [], [] # true
        _energy, _force, _stress = [], [], [] # predicted
 
        for i in range(len(list(db.keys()))):
            no_of_atoms = len(db[str(i)]['force'])
            Energy, Force, Stress = self.calculate_properties(db[str(i)],       # Energy per atom
                                    self.force_coefficient, self.stress_coefficient)
            
            # Store energy into list
            true_energy = db[str(i)]['energy'] / no_of_atoms
            energy.append(true_energy)
            _energy.append(Energy.sum() / no_of_atoms)

            if self.force_coefficient:
                true_force = np.ravel(db[str(i)]['force'])
                Force = np.ravel(Force)
                for m in range(len(true_force)):
                    force.append(true_force[m])
                    _force.append(Force[m])

            if self.stress_coefficient and (db[str(i)]['group'] in self.stress_group):
                true_stress = np.ravel(db[str(i)]['stress'])#.flat[[0,3,5,3,1,4,5,4,2]]
                Stress = np.ravel(Stress)
                for m in range(len(true_stress)):
                    stress.append(true_stress[m])
                    _stress.append(Stress[m])

        energy, force, stress = np.asarray(energy), np.asarray(force), np.asarray(stress)
        _energy, _force, _stress = np.asarray(_energy), np.asarray(_force), np.asarray(_stress)

        # Dump the true and predicted values into text file.
        self.dump_evaluate(_energy, energy, filename=figname[:-4]+'Energy.txt')
        if self.force_coefficient:
            self.dump_evaluate(_force, force, filename=figname[:-4]+'Force.txt')
        if self.stress_coefficient:
            self.dump_evaluate(_stress, stress, filename=figname[:-4]+'Stress.txt')
        
        # Calculate the statistical metrics for energy.
        E_mae = self.mean_absolute_error(energy, _energy)
        E_mse = self.mean_squared_error(energy, _energy)
        E_r2 = self.r2_score(energy, _energy)
        print("The results for energy: ")
        print("    Energy R2     {:8.6f}".format(E_r2))
        print("    Energy MAE    {:8.6f}".format(E_mae))
        print("    Energy RMSE   {:8.6f}".format(E_mse))

        # Plotting the energy results.
        energy_str = 'Energy: r2({:.4f}), MAE({:.4f} {}/atom)'. \
                     format(E_r2, E_mae, self.unit)
        plt.title(energy_str)
        plt.scatter(energy, _energy, label='Energy', s=5)
        plt.legend(loc=2)
        plt.xlabel('True ({}/atom)'.format(self.unit))
        plt.ylabel('Prediction ({}/atom)'.format(self.unit))
        plt.tight_layout()
        plt.savefig(self.path+'Energy_'+figname)
        plt.close()
        print("The energy figure is exported to: {:s}".format(self.path+'Energy_'+figname))
        print("\n")

        if self.force_coefficient:
            F_mae = self.mean_absolute_error(force, _force)
            F_mse = self.mean_squared_error(force, _force)
            F_r2 = self.r2_score(force, _force)
            print("The results for force: ")
            print("    Force R2      {:8.6f}".format(F_r2))
            print("    Force MAE     {:8.6f}".format(F_mae))
            print("    Force RMSE    {:8.6f}".format(F_mse))

            # Plotting the forces results.
            length = 'A'
            if self.unit == 'Ha':
                length == 'Bohr'
            force_str = 'Force: r2({:.4f}), MAE({:.3f} {}/{})'. \
                        format(F_r2, F_mae, self.unit, length)
            plt.title(force_str)
            plt.scatter(force, _force, s=5, label='Force')
            plt.legend(loc=2)
            plt.xlabel('True ({}/{})'.format(self.unit, length))
            plt.ylabel('Prediction ({}/{})'.format(self.unit, length))
            plt.tight_layout()
            plt.savefig(self.path+'Force_'+figname)
            plt.close()
            print("The force figure is exported to: {:s}".format(self.path+'Force_'+figname))
            print("\n")

        else:
            F_mae, F_mse, F_r2 = None, None, None

        if self.stress_coefficient:
            S_mae = self.mean_absolute_error(stress, _stress)
            S_mse = self.mean_squared_error(stress, _stress)
            S_r2 = self.r2_score(stress, _stress)
            print("The results for stress: ")
            print("    Stress R2      {:8.6f}".format(S_r2))
            print("    Stress MAE     {:8.6f}".format(S_mae))
            print("    Stress RMSE    {:8.6f}".format(S_mse))

            # Plotting the stress results.
            stress_str = 'Stress: r2({:.4f}), MAE({:.3f} GPa)'. \
                        format(S_r2, S_mae)
            plt.title(stress_str)
            plt.scatter(stress, _stress, s=5, label='Stress')
            plt.legend(loc=2)
            plt.xlabel('True (GPa)')
            plt.ylabel('Prediction (GPa)')
            plt.tight_layout()
            plt.savefig(self.path+'Stress_'+figname)
            plt.close()
            print("The stress figure is exported to: {:s}".format(self.path+'Stress_'+figname))
            print("\n")
        else:
            S_mae, S_mse, S_r2 = None, None, None

        return (E_mae, E_mse, E_r2, F_mae, F_mse, F_r2, S_mae, S_mse, S_r2)


    def LinearRegression(self, X, y, w, alpha, norm=2):
        """ Perform linear regression. """
        
        m = X.shape[1] # The shape of the descriptors

        _X = X * np.sqrt(np.expand_dims(w, axis=1))
        _y = y * np.sqrt(w)

        if self.alpha:
            if norm == 1:
                theta = np.linalg.lstsq(_X.T.dot(_X), _X.T.dot(_y) - alpha, rcond=None)[0]
            elif norm == 2:
                theta = np.linalg.lstsq(_X.T.dot(_X) + alpha * np.identity(m), _X.T.dot(_y), rcond=None)[0]
            else:
                msg = f"Regularization with {norm} norm is not implemented yet."
                raise NotImplementedError(msg)
        else:
            theta = np.linalg.lstsq(_X, _y, rcond=None)[0]

        return theta


    def save_checkpoint(self, des_info, filename=None):
        """ Save Polynomial Regression model to PyTorch. """
        _filename = self.path

        if filename:
            _filename += filename
        else:
            _filename += 'PolyReg-checkpoint.pth'

        checkpoint = {'elements': self.elements,
                      'algorithm': 'PR',
                      'force_coefficient': self.force_coefficient,
                      'path': self.path,
                      'quadratic': self.quadratic,
                      'coef_': self.coef_,
                      'des_info': des_info}

        save(checkpoint, _filename)
        if des_info['type'] in ['SNAP', 'snap', 'SO3', 'SOAP']:
            self.save_weights_to_txt(des_info)
        print("The Linear Regression Potential is exported to {:s}".format(_filename))
        print("\n")


    def save_weights_to_txt(self, des_info):
        """ Saving the model weights to txt file. """
        with open(self.path+"PR_weights.txt", "w") as f:
            f.write("# Polynomial Regression weights generated in PyXtal_FF \n")
            f.write("# total_species ncoefficient \n\n")
            f.write(f"{len(self.elements)} {self.d_max+1} \n")
            count = 0
            for element in self.elements:
                #if des_info['type'] in ['SNAP', 'snap']:
                #    f.write(f"{element} 0.5 {des_info['weights'][element]} \n")
                #else:
                #    f.write(f"{element} \n")
                for _ in range(self.d_max+1):
                    f.write(f"{self.coef_[count]} \n")
                    count += 1

        with open(self.path+"DescriptorParam.txt", "w") as f:
            f.write("# Descriptor parameters generated in PyXtal_FF \n\n")
            f.write("# Required \n")
            f.write(f"rcutfac {des_info['Rc']} \n")
            
            if des_info['type'] in ['SO3', 'SOAP']:
                f.write(f"nmax {des_info['parameters']['nmax']} \n")
                f.write(f"lmax {des_info['parameters']['lmax']} \n")
                f.write(f"alpha {des_info['parameters']['alpha']} \n\n")
            else:
                f.write(f"twojmax {des_info['parameters']['lmax']*2} \n\n")

            f.write("# Elements \n\n")
            f.write(f"nelems {len(self.elements)} \n")
            f.write("elems ")
            for element in self.elements:
                f.write(f"{element} ")
            f.write("\n")

            if des_info['type'] in ['snap', 'SNAP', 'SO3', 'SOAP']:
                f.write("radelems ")
                for element in self.elements:
                    f.write("0.5 ")
                f.write("\n")

                if des_info['type'] in ['snap', 'SNAP']:
                    f.write("welems ")
                    for element in self.elements:
                        f.write(f"{des_info['weights'][element]} ")
                    f.write("\n\n")
                else:
                    f.write("welems ")
                    ele = Element(self.elements)
                    atomic_numbers = ele.get_Z()
                    for num in atomic_numbers:
                        f.write(f"{num} ")
                    f.write("\n\n")

            if des_info['type'] in ['snap', 'SNAP']:
                f.write(f"rfac0 {des_info['parameters']['rfac']} \n")
                f.write(f"rmin0 0 ")
                f.write("\n")
                f.write("switchflag 1 \n")
                f.write("bzeroflag 0 \n")


    def load_checkpoint(self, filename=None):
        """ Load Polynomial Regression file from PyTorch. """
        checkpoint = load(filename, weights_only=False) # New in PyTorch 2.6

        # Inconsistent algorithm.
        if checkpoint['algorithm'] != 'PR':
            msg = "The loaded algorithm is not Polynomial Regression."
            raise NotImplementedError(msg)
        
        # Check the consistency with the system of elements
        msg = f"The system, {self.elements}, are not consistent with "\
                    +"the loaded system, {checkpoint['elements']}."

        if len(self.elements) != len(checkpoint['elements']):
            raise ValueError(msg)
        
        for i in range(len(self.elements)):
            if self.elements[i] != checkpoint['elements'][i]:
                raise ValueError(msg)
        
        self.coef_ = checkpoint['coef_']
        self.quadratic = checkpoint['quadratic']

        return checkpoint['des_info']


    def calculate_properties(self, descriptor, bforce=True, bstress=False):
        """ A routine to compute energy, forces, and stress.
        
        Parameters:
        -----------
        descriptor: list
            list of x, dxdr, and rdxdr.
        benergy, bforce, bstress: bool
            If False, excluding the property from calculation.

        Returns:
        --------
        energy: float
            The predicted energy
        forces: 2D array [N_atom, 3] (if dxdr is provided)
            The predicted forces
        stress: 2D array [3, 3] (if rdxdr is provided)
            The predicted stress
        """
        no_of_atoms = len(descriptor['elements'])
        energies, force, stress = np.zeros([no_of_atoms]), np.zeros([no_of_atoms, 3]), np.zeros([6])
        
        X = self.parse_descriptors({'0': descriptor}, fc=bforce, sc=bstress, train=False)
        
        _y = np.dot(X, self.coef_) # Calculate properties

        energies = _y[:no_of_atoms]

        #energy = _y[0] / no_of_atoms # get energy/atom
        
        if bforce: # get force
            force += np.reshape(_y[no_of_atoms:no_of_atoms+(no_of_atoms*3)], (no_of_atoms, 3))

        if bstress: # get stress
            stress += _y[-6:]*eV2GPa # in GPa
        
        return energies, force, stress


    def mean_absolute_error(self, true, predicted):
        """ Calculate mean absolute error of energy or force. """
        return sum(abs(true-predicted)/len(true))


    def mean_squared_error(self, true, predicted):
        """ Calculate mean square error of energy or force. """
        return np.sqrt(sum((true-predicted) ** 2 /len(true)))


    def r2_score(self, true, predicted):
        """ Calculate the r square of energy or force. """
        t_bar = sum(true)/len(true)
        square_error = sum((true-predicted) ** 2)
        true_variance = sum((true-t_bar) ** 2)
        return 1 - square_error / true_variance


    def dump_evaluate(self, predicted, true, filename):
        """ Dump the evaluate results to text files. """
        absolute_diff = np.abs(np.subtract(predicted, true))
        combine = np.vstack((predicted, true, absolute_diff)).T
        np.savetxt(self.path+filename, combine, header='Predicted True Diff', fmt='%.7e')


    def parse_descriptors(self, data, fc=True, sc=False, train=True):
        """ Parse descriptors and its gradient to 2-D array. 
        
        Returns
        -------
        X: 2-D array [n+m*3, d]
            d is the total number of descriptors, n is the total
            number of structures, and m is the total number atoms
            in the entire structures. If force_coefficient is None,
            X has the shape of [n, d].
        """
        if train:
            db = shelve.open(self.path+data)
            no_of_structures = self.no_of_structures
            no_of_atoms = self.no_of_atoms
            stress_components = self.stress_components
        else:
            db = data
            no_of_structures = 1 # 1 for train is false
            no_of_atoms = len(data['0']['elements']) if fc else 0
            stress_components = 6 if sc else 0
        
        # Determine the total number of descriptors based on SNAP or qSNAP.
        # Note: d_max != self.d_max
        if self.d_max is None: #enable calculator works
            self.d_max = len(db['0']['x'][0])

        if self.quadratic:
            d_max = (self.d_max**2+3*self.d_max)//2 # (d^2 + 3*d) / 2
        else:
            d_max = self.d_max

        columns = (1+d_max)*len(self.elements)
        rows = no_of_structures if train else no_of_atoms
        rows += no_of_atoms * 3 if fc else 0 # x, y, and z
        rows += stress_components if sc else 0 # xx, xy, xz, ..., zz
        
        X = np.zeros([rows, columns])

        # Fill in X.
        xcount = 0
        for i in range(no_of_structures):
            data = db[str(i)]
            if train:
                try:
                    _group = True if data['group'] in self.stress_group else False
                except:
                    _group = False
            else:
                _group = True
           
            if self.quadratic:
                _x = data['x'][:, :self.d_max]
                x = np.zeros((len(_x), d_max))
                x[:, :self.d_max] += data['x'][:, :self.d_max]
                if fc:
                    _dxdr = data['dxdr'][:, :self.d_max, :]
                    dxdr = np.zeros([_dxdr.shape[0], d_max, 3])
                    dxdr[:, :self.d_max, :] += _dxdr
                    seq = data['seq']
                
                #if train:
                if sc and _group: #(data['group'] in self.stress_group):
                    _rdxdr = data['rdxdr'][:, :self.d_max, :]
                    rdxdr = np.zeros([len(_x), d_max, 6])
                    rdxdr[:, :self.d_max, :] += _rdxdr
                
                # self-term for x, dxdr, and rdxdr
                x_square = 0.5 * _x ** 2
                if fc:
                    dxdr_square = np.zeros_like(_dxdr)
                    for i in range(len(_x)):
                        arr = np.where(seq[:, 0]==i)[0]
                        dxdr_square[arr] = np.einsum('jkl, k->jkl', _dxdr[arr], _x[i])
                    #dxdr_square = np.einsum('ijkl,ik->ijkl', _dxdr, _x)
                if sc and _group:
                    rdxdr_square = np.einsum('ijk,ij->ijk', _rdxdr, _x)
                
                dcount = self.d_max
                for d1 in range(self.d_max):
                    # Cross term for x, dxdr, and rdxdr
                    x_d1_d2 = np.einsum('i, ij->ij', _x[:, d1], _x[:, d1+1:])
                    if fc:
                        shp = _dxdr.shape
                        dxdr_d1_d2 = np.zeros([shp[0], shp[1]-1-d1, shp[2]])
                        for i in range(len(_x)):
                            arr = np.where(seq[:, 0]==i)[0]
                            dxdr_d1_d2[arr] = np.einsum('ij,k->ikj', _dxdr[arr][:, d1, :], _x[i, d1+1:]) + \
                                    (_dxdr[arr][:, d1+1:, :] * _x[i, d1])

                        #dxdr_d1_d2 = np.einsum('ijl,ik->ijkl', _dxdr[:, :, d1, :], _x[:, d1+1:]) + \
                                #             np.einsum('ijkl,i->ijkl', _dxdr[:, :, d1+1:, :], _x[:, d1])
                    if sc and _group:
                        rdxdr_d1_d2 = np.einsum('ik,ij->ijk', _rdxdr[:, d1, :], _x[:, d1+1:]) + \
                                      np.einsum('ijk,i->ijk', _rdxdr[:, d1+1:, :], _x[:, d1])
                    
                    # Append for x, dxdr, rdxdr
                    x[:, dcount] += x_square[:, d1]
                    if fc:
                        dxdr[:, dcount, :] += dxdr_square[:, d1, :]
                    if sc and _group:
                        rdxdr[:, dcount, :] += rdxdr_square[:, d1, :]

                    dcount += 1
                    
                    x[:, dcount:dcount+len(x_d1_d2[0])] += x_d1_d2
                    if fc:
                        dxdr[:, dcount:dcount+len(x_d1_d2[0]), :] += dxdr_d1_d2
                    if sc and _group:
                        rdxdr[:, dcount:dcount+len(x_d1_d2[0]), :] += rdxdr_d1_d2
                    dcount += len(x_d1_d2[0])
                
            else:
                x = data['x'][:, :d_max]
                if fc:
                    seq = data['seq']   
                    dxdr = data['dxdr'][:, :d_max, :]
                if sc and _group:
                    rdxdr = data['rdxdr'][:, :d_max, :]
            
            elements = data['elements']
            
            # Arranging x and dxdr for energy and forces.
            bias_weights = 1.0
            
            if train:
                sna = np.zeros([len(self.elements), 1+d_max])
                if fc:
                    snad = np.zeros([len(self.elements), len(x), 1+d_max, 3])
                if sc and _group:
                    snav = np.zeros([len(self.elements), 1+d_max, 6])
                
                _sna, _snad, _snav, _count = {}, {}, {}, {}
                for element in self.elements:
                    _sna[element] = None
                    _snad[element] = None
                    _snav[element] = None
                    _count[element] = 0
                
                # Loop over the number of atoms in a structure.
                for e, element in enumerate(elements):
                    if _sna[element] is None:
                        _sna[element] = 1 * x[e]
                        if fc:
                            shp = snad.shape
                            _snad[element] = np.zeros([shp[1], shp[2]-1, shp[3]])
                            arr = np.where(seq[:, 0]==e)[0]
                            _seq = seq[arr][:, 1]
                            _snad[element][_seq] = -1 * dxdr[arr]
                        if sc and _group:
                            _snav[element] = -1 * rdxdr[e]  # [d, 6]
                    else:
                        _sna[element] += x[e]
                        if fc:
                            arr = np.where(seq[:, 0]==e)[0]
                            _seq = seq[arr][:, 1]
                            _snad[element][_seq] -= dxdr[arr]
                        if sc and _group: 
                            _snav[element] -= rdxdr[e]
                    _count[element] += 1

                for e, element in enumerate(self.elements):
                    if _count[element] > 0:
                        #_sna[element] /= _count[element]
                        sna[e, :] += np.hstack(([bias_weights*_count[element]], _sna[element]))
                        if fc:
                            snad[e, :, 1:, :] += _snad[element]
                        if sc and _group:
                            snav[e, 1:, :] += _snav[element]
                        
                # X for energy
                X[xcount, :] += sna.ravel()
                xcount += 1

                # X for forces.
                if fc:
                    for j in range(snad.shape[1]):
                        for k in range(snad.shape[3]):
                            X[xcount, :] += snad[:, j, :, k].ravel()
                            xcount += 1
                
                # X for stress.
                if sc and _group: 
                    shp = snav.shape
                    X[xcount:xcount+6, :] = snav.reshape([shp[0]*shp[1], shp[2]]).T
                    xcount += 6

            else:
                if fc:
                    snad = np.zeros([len(self.elements), len(x), 1+d_max, 3])
                if sc and _group:
                    snav = np.zeros([len(self.elements), 1+d_max, 6])
                _snad, _snav, _count = {}, {}, {}
                for element in self.elements:
                    _snad[element] = None
                    _snav[element] = None
                    _count[element] = 0

                # Loop over the number of atoms in a structure.
                for e, element in enumerate(elements):
                    elem_cnt = self.elements.index(element)
                    X[xcount, elem_cnt*(1+d_max):(elem_cnt+1)*(1+d_max)] = np.hstack(([bias_weights], x[e]))
                    xcount += 1
                    if _snad[element] is None:
                        if fc:
                            shp = snad.shape
                            _snad[element] = np.zeros([shp[1], shp[2]-1, shp[3]])
                            arr = np.where(seq[:, 0]==e)[0]
                            _seq = seq[arr][:, 1]
                            _snad[element][_seq] = -1 * dxdr[arr]
                        if sc and _group:
                            _snav[element] = -1 * rdxdr[e]  # [d, 6]
                    else:
                        if fc:
                            arr = np.where(seq[:, 0]==e)[0]
                            _seq = seq[arr][:, 1]
                            _snad[element][_seq] -= dxdr[arr]
                        if sc and _group: 
                            _snav[element] -= rdxdr[e]
                    _count[element] += 1

                for e, element in enumerate(self.elements):
                    if _count[element] > 0:
                        if fc:
                            snad[e, :, 1:, :] += _snad[element]
                        if sc and _group:
                            snav[e, 1:, :] += _snav[element]

                # X for forces.
                if fc:
                    for j in range(snad.shape[1]):
                        for k in range(snad.shape[3]):
                            X[xcount, :] += snad[:, j, :, k].ravel()
                            xcount += 1
                
                # X for stress.
                if sc and _group: 
                    shp = snav.shape
                    X[xcount:xcount+6, :] = snav.reshape([shp[0]*shp[1], shp[2]]).T
                    xcount += 6
        
        if train:
            db.close()

        return X


    def parse_features(self, data):
        """ Parse features (energy, forces, and stress) into 1-D array.
        
        Returns
        -------
        y: 1-D array [n+m*3+n*3*3,]
            y contains the energy, forces, and stress of structures 
            in 1-D array. Force and stress may not be present.
            
            n = # of structures
            m = # of atoms in a unit cell

        w: 1-D array [n+m*3+n*3*3,]
            w contains the relative importance between energy, forces, 
            and stress.
        """
        db = shelve.open(self.path+data)
        self.no_of_atoms = 0
        self.stress_components = 0

        y = None # store the features (energy+forces+stress)
        w = None # weight of each sample
        
        for i in range(self.no_of_structures):
            data = db[str(i)]
            energy = np.array([data['energy']])
            w_energy = np.array([1.])

            if self.force_coefficient:
                force = np.array(data['force']).ravel()
                w_force = np.array([self.force_coefficient]*len(force))

                if self.stress_coefficient and (data['group'] in self.stress_group):     # energy + forces + stress
                    stress = np.array(data['stress']).ravel()#.flat[[0,3,5,3,1,4,5,4,2]]
                    w_stress = np.array([self.stress_coefficient]*len(stress))
                    self.stress_components += 6
                    
                    if y is None:
                        y = np.concatenate((energy, force, stress))
                        w = np.concatenate((w_energy, w_force, w_stress))
                    else:
                        y = np.concatenate((y, energy, force, stress))
                        w = np.concatenate((w, w_energy, w_force, w_stress))
                else:                                                                           # energy + forces
                    if y is None:
                        y = np.concatenate((energy, force))
                        w = np.concatenate((w_energy, w_force))
                    else:
                        y = np.concatenate((y, energy, force))
                        w = np.concatenate((w, w_energy, w_force))
                
                # Count the number of atoms for the entire structures.
                self.no_of_atoms += len(data['force'])

            else:
                if self.stress_coefficient and (data['group'] in self.stress_group):    # energy + stress
                    stress = np.array(data['stress'])#.flat[[0,3,5,3,1,4,5,4,2]]
                    w_stress = np.array([self.stress_coefficient]*len(stress))
                    self.stress_components += 6
                    
                    if y is None:
                        y = np.concatenate((energy, stress))
                        w = np.concatenate((w_energy, w_stress))
                    else:
                        y = np.concatenate((y, energy, stress))
                        w = np.concatenate((w, w_energy, w_stress))

                else:                                                                           # energy only
                    if y is None:
                        y = energy
                        w = w_energy
                    else:
                        y = np.concatenate((y, energy))
                        w = np.concatenate((w, w_energy))
        
        db.close()
        gc.collect()

        return y, w




class DescriptorSO3:
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

    def __init__(self, nmax=3, lmax=3, rcut=3.5, alpha=2.0, derivative=True, stress=False, cutoff_function='cosine', weight_on=False, primitive=False):
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

    def __repr__(self):
        return str(self)

    def load_from_dict(self, dict0):
        self.nmax = dict0["nmax"]
        self.lmax = dict0["lmax"]
        self.rcut = dict0["rcut"]
        self.alpha = dict0["alpha"]
        self.derivative = dict0["derivative"]
        self.stress = dict0["stress"]

    def save_dict(self):
        """
        save the model as a dictionary in json
        """
        dict = {"nmax": self.nmax,
                "lmax": self.lmax,
                "rcut": self.rcut,
                "alpha": self.alpha,
                "derivative": self.derivative,
                "stress": self.stress,
                "_type": "SO3",
               }
        return dict

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
    def derivative(self):
        return self._derivative

    @derivative.setter
    def derivative(self, derivative):
        if isinstance(derivative, bool) is True:
            self._derivative = derivative
        else:
            raise ValueError('derivative must be a boolean value')

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
    def cutoff_function(self):
        return self._cutoff_function

    @cutoff_function.setter
    def cutoff_function(self, cutoff_function):
        if isinstance(cutoff_function, str) is True:
            # more conditions
            if cutoff_function == 'cosine':
                self._cutoff_function = Cosine
            elif cutoff_function == 'tanh':
                self._cutoff_function = Tanh
            elif cutoff_function == 'poly1':
                self._cutoff_function = Poly1
            elif cutoff_function == 'poly2':
                self._cutoff_function = Poly2
            elif cutoff_function == 'poly3':
                self._cutoff_function = Poly3
            elif cutoff_function == 'poly4':
                self._cutoff_function = Poly4
            elif cutoff_function == 'exp':
                self._cutoff_function = Exponent
            elif cutoff_function == 'unity':
                self._cutoff_function = Unity
            else:
                raise NotImplementedError('The requested cutoff function has not been implemented')
        else:
            raise ValueError('You must specify the cutoff function as a string')

    def clear_memory(self):
        '''
        Clears all memory that isn't an essential attribute for the calculator
        '''
        attrs = list(vars(self).keys())
        for attr in attrs:
            if attr not in {'_nmax', '_lmax', '_rcut', '_alpha', '_derivative', '_stress', '_cutoff_function', 'weight_on', 'primitive'}:
                delattr(self, attr)
        return

    def calculate(self, atoms, atom_ids=None):
        '''
        Calculates the SO(3) power spectrum components of the
        smoothened atomic neighbor density function
        for given nmax, lmax, rcut, and alpha.

        Args:
            atoms: an ASE atoms object corresponding to the desired
                   atomic arrangement

            backend: string, specifies the method to compute the neighborlist
                     elements, either ASE or pymatgen
        '''
        self._atoms = atoms

        self.build_neighbor_list(atom_ids)
        self.initialize_arrays()

        ncoefs = self.nmax*(self.nmax+1)//2*(self.lmax+1)
        tril_indices = np.tril_indices(self.nmax, k=0)

        ls = np.arange(self.lmax+1)
        norm = np.sqrt(2*np.sqrt(2)*np.pi/np.sqrt(2*ls+1))

        if self.derivative:
            # get expansion coefficients and derivatives
            cs, dcs = compute_dcs(self.neighborlist, self.nmax, self.lmax, self.rcut, self.alpha, self._cutoff_function)
            # weight cs and dcs
            cs *= self.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis]
            dcs *= self.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis,np.newaxis]
            cs = np.einsum('inlm,l->inlm', cs, norm)
            dcs = np.einsum('inlmj,l->inlmj', dcs, norm)
            Ris = self.center_atoms
            Rjs = self.neighborlist + Ris
            for i in np.unique(self.seq[:,0]):
                # find atoms for which i is the center
                centers = self.neighbor_indices[:,0] == i
                # find neighbors for which i is not the index
                neighs = self.neighbor_indices[:,1] != i
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
                ijs = self.neighbor_indices[centers]
                # loop over unique neighbor indices
                for j in np.unique(ijs[:,1]):
                    # get the location of ij pairs in the NL
                    # and therefore dP
                    ijlocs = self.neighbor_indices[centers,1] == j
                    # get the location of the dplist element
                    temp = self.seq == np.array([i,j])
                    seqloc = temp[:,0]*temp[:,1]
                    # sum over ij pairs
                    dPsum = np.sum(dP[ijlocs], axis=0)
                    rdPjsum = np.sum(rdPj[ijlocs], axis=0)
                    # flatten into dplist and rdplist
                    self._dplist[seqloc] += (dPsum[tril_indices].flatten()).reshape(ncoefs,3)
                    self._pstress[seqloc] -= (rdPjsum[tril_indices].flatten()).reshape(ncoefs,3,3)

                # get unique elements and store in feature vector
                self._plist[i] = P[tril_indices].flatten()
                # get location if ii pair in seq
                temp = self.seq == np.array([i,i])
                iiloc = temp[:,0]*temp[:,1]
                # get location of all ijs in seq
                ilocs = self.seq[:,0] == i
                self._dplist[iiloc] -= np.sum(self._dplist[ilocs],axis=0)
                rdPisum = np.sum(rdPi, axis=0)
                self._pstress[iiloc] += (rdPisum[tril_indices].flatten()).reshape(ncoefs,3,3)


            x = {'x':self._plist, 'dxdr':self._dplist,
                 'elements':list(atoms.symbols), 'seq':self.seq}
            if self._stress:
                vol = atoms.get_volume()
                x['rdxdr'] = -self._pstress/vol
            else:
                x['rdxdr'] = None

        else:
            cs = compute_cs(self.neighborlist, self.nmax, self.lmax, self.rcut, self.alpha, self._cutoff_function)
            cs *= self.atomic_weights[:,np.newaxis,np.newaxis,np.newaxis]
            cs = np.einsum('inlm,l->inlm', cs, norm)
            # everything good up to here
            for i in np.unique(self.seq[:,0]):
                centers = self.neighbor_indices[:,0] == i
                ctot = cs[centers].sum(axis=0)
                P = np.einsum('ijk,ljk->ilj', ctot, np.conj(ctot)).real
                self._plist[i] = P[tril_indices].flatten()
            x = {'x':self._plist, 'dxdr': None, 'rdxdr': None, 'elements':list(atoms.symbols)}

        self.clear_memory()
        return x

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
        self._dplist = np.zeros((len(self.seq), ncoefs, 3), dtype=np.float64)
        self._pstress = np.zeros((len(self.seq), ncoefs, 3, 3), dtype=np.float64)

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
            nl = PrimitiveNeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
            nl.build(atoms.pbc, atoms.cell, atoms.get_scaled_positions())
        else:
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

def Cosine(Rij, Rc, derivative=False):
    # Rij is the norm
    if derivative is False:
        result = 0.5 * (np.cos(np.pi * Rij / Rc) + 1.)
    else:
        result = -0.5 * np.pi / Rc * np.sin(np.pi * Rij / Rc)
    return result

def Tanh(Rij, Rc, derivative=False):

    if derivative is False:
        result = np.tanh(1-Rij/Rc)**3

    else:
        tanh_square = np.tanh(1-Rij/Rc)**2
        result = - (3/Rc) * tanh_square * (1-tanh_square)
    return result

def Poly1(Rij, Rc, derivative=False):

    if derivative is False:
        x = Rij/Rc
        x_square = x**2
        result = x_square * (2*x-3) + 1

    else:
        term1 = (6 / Rc**2) * Rij
        term2 = Rij/Rc - 1
        result = term1*term2
    return result

def Poly2(Rij, Rc, derivative=False):

    if derivative is False:
        x = Rij/Rc
        result = x**3 * (x*(15-6*x)-10) + 1

    else:
        x = Rij/Rc
        result = (-30/Rc) * (x**2 * (x-1)**2)
    return result

def Poly3(Rij, Rc, derivative=False):

    if derivative is False:
        x = Rij/Rc
        result = x**4*(x*(x*(20*x-70)+84)-35)+1

    else:
        x = Rij/Rc
        result = (140/Rc) * (x**3 * (x-1)**3)
    return result

def Poly4(Rij, Rc, derivative=False):

    if derivative is False:
        x = Rij/Rc
        result = x**5*(x*(x*(x*(315-70*x)-540)+420)-126)+1

    else:
        x = Rij/Rc
        result = (-630/Rc) * (x**4 * (x-1)**4)
    return result

def Exponent(Rij, Rc, derivative=False):

    if derivative is False:
        x = Rij/Rc
        try:
            result = np.exp(1 - 1/(1-x**2))
        except:
            result = 0

    else:
        x = Rij/Rc
        try:
            result = 2*x * np.exp(1 - 1/(1-x**2)) / (1+x**2)**2
        except:
            result = 0
            return result

def Unity(Rij, Rc, derivative=False):
    if derivative is False:
        return np.ones(len(Rij))

    else:
        return np.ones(len(Rij))

def W(nmax):
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

def phi(r, alpha, rcut):
    '''
    See g below
    '''
    return (rcut-r)**(alpha+2)/np.sqrt(2*rcut**(2*alpha+7)/(2*alpha+5)/(2*alpha+6)/(2*alpha+7))

def g(r, n, nmax, rcut, w):

    Sum = 0.0
    for alpha in range(1, nmax+1):
        Sum += w[n-1, alpha-1]*phi(r, alpha, rcut)

    return Sum

def GaussChebyshevQuadrature(nmax,lmax):
    NQuad = (nmax+lmax+1)*10
    quad_array = np.zeros(NQuad, dtype=np.float64)
    for i in range(1,NQuad+1,1):
        # roots of Chebyshev polynomial of degree N
        x = np.cos((2*i-1)*np.pi/2/NQuad)
        quad_array[i-1] = x
    return quad_array, np.pi/NQuad

def compute_cs(pos, nmax, lmax, rcut, alpha, cutoff):

    # compute the overlap matrix
    w = W(nmax)

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
            ylms[:,l,midx] = sph_harm(m, l, phis, thetas)

    # multiply the spherical harmonics and the radial inner product
    Y_mul_innerprod = np.einsum('ijk,ilj->iljk', ylms, integral_array)

    # multiply the gaussians into the expression
    C = np.einsum('i,ijkl->ijkl', exparray, Y_mul_innerprod)
    return C

def compute_dcs(pos, nmax, lmax, rcut, alpha, cutoff):
    # compute the overlap matrix
    w = W(nmax)

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
            ylms[:,l,midx] = sph_harm(m, l, phis, thetas)


    for l in range(1, lmax+1):
        for m in range(-l, l+1, 1):
            midx = msize//2 + m
            Midx = Msize//2 + m
            # get gradient with recpect to spherical covariant components
            xcov0 = -np.sqrt(((l+1)**2-m**2)/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx]/Ris

            if abs(m) <= l-1:
                xcov0 += np.sqrt((l**2-m**2)/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx]/Ris


            xcovpl1 = -np.sqrt((l+m+1)*(l+m+2)/2/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx+1]/Ris

            if abs(m+1) <= l-1:
                xcovpl1 -= np.sqrt((l-m-1)*(l-m)/2/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx+1]/Ris


            xcovm1 = -np.sqrt((l-m+1)*(l-m+2)/2/(2*l+1)/(2*l+3))*l*ylms[:,l+1,Midx-1]/Ris

            if abs(m-1) <= l-1:
                xcovm1 -= np.sqrt((l+m-1)*(l+m)/2/(2*l-1)/(2*l+1))*(l+1)*ylms[:,l-1,Midx-1]/Ris

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







class Database():
    def __init__(self, name):
        self.name = name
        print("name = ", name)
        self.database = shelve.open(self.name)
        
        self.msg = "Index must be an integer in the interval [0,len(self)]"
        self.length = len(list(self.database.keys()))


    def __len__(self):
        return self.length


    def __setitem__(self, index, value):
        if isinstance(index, int) and index >= 0:
            self.database[str(index)] = value
        else:
            raise IndexError(self.msg)


    def __getitem__(self, index):
        if isinstance(index, int) and index >= 0 and index < len(self):
            return self.database[str(index)]
        else:
            raise IndexError(self.msg)


    def __delitem__(self, index):
        if isinstance(index, int) and index >= 0:
            del(self.database[str(index)])
        else:
            raise IndexError(self.msg)


    def insert(self, index, value):
        """ Insert the value to the dictionary at index. """
        if isinstance(index, int) and index >= 0:
            pass
        else:
            raise IndexError(self.msg)
        self[index] = value


    def append(self, value):
        """ Append value to the end of the sequence. """
        self.insert(len(self), value)


    def close(self):
        self.database.close()


    def store(self, structure_file, function, storage, ase_db=None):
        """ Map structures to descriptors and store them, including features, to database.
        If compute is False, print pre-computed descriptors message. """
        if function['base_potential']:
            raise NotImplementedError("base_potential is not implemented")
            #self.base_potential = ZBL(function['base_potential']['inner'],
            #                          function['base_potential']['outer'])
        else:
            self.base_potential = None

        if storage:
            if os.path.isdir(structure_file):
                fmt = 'dat'
            elif structure_file.find('json') > 0:
                fmt = 'json'
            elif structure_file.find('xyz') > 0:
                fmt = 'xyz'
            elif structure_file.find('db') > 0:
                fmt = 'db'
            elif structure_file.find('traj') > 0:
                fmt = 'traj'
            else:
                fmt = 'vasp-out'
        
            # extract the structures and energy, forces, and stress information.
            if os.path.exists(structure_file):
                if fmt == 'xyz':
                    data = parse_xyz(structure_file)
                elif fmt == 'db':
                    data = parse_ase_db(structure_file)
                else:
                    raise NotImplementedError
            else:
                raise FileNotFoundError(structure_file + ' cannot be found from the given path.')
            print("{:d} structures have been loaded.".format(len(data)))
            
            self.add(function, data)

            if ase_db is not None and fmt != 'db':
                convert_to_ase_db(data, ase_db)
                print("save the structures to {:}.".format(ase_db))

        else:
            print(f"Features and precomputed descriptors exist: {self.name}.dat\n")


    def add(self, function, data):
        """ Add descriptors for all structures to database. """
        print('Computing the descriptors...')

        _N = deepcopy(function['N'])
        _cpu = 1

        N1 = len(data)
        if _N is not None and _N < N1:
            lists = range(_N)
        else:
            lists = range(N1)
        
        for i, index in enumerate(lists):
            d = self.compute(function, data[index])
            self.append(d)
            self.length += 1
            print('\r{:4d} out of {:4d}'.format(i+1, len(lists)), flush=True, end='')

        print(f"\nSaving descriptor-feature data to {self.name}.dat\n")


    def compute(self, function, data):
        """ Compute descriptor for one structure to the database. """

        if function['type'] in ['BehlerParrinello', 'ACSF']:
            raise NotImplementedError

        elif function['type'] in ['wACSF', 'wacsf']:
            raise NotImplementedError
        
        elif function['type'] in ['SO4', 'Bispectrum', 'bispectrum']:
            raise NotImplementedError
        
        elif function['type'] in ['SO3', 'SOAP', 'soap']:
            d = DescriptorSO3(
                    function['parameters']['nmax'],
                    function['parameters']['lmax'],
                    function['Rc'],
                    alpha=function['parameters']['alpha'],
                    derivative=function['force'],
                    stress=function['stress']
                ).calculate(data['structure'])

        elif function['type'] in ['EAD', 'ead']:
            raise NotImplementedError
        
        elif function['type'] in ['SNAP', 'snap']:
            raise NotImplementedError

        else:
            msg = f"{function['type']} is not implemented"
            raise NotImplementedError(msg)

        if d['rdxdr'] is not None:
            N = d['x'].shape[0]
            L = d['x'].shape[1]
            rdxdr = np.zeros([N, L, 3, 3])
            for _m in range(N):
                ids = np.where(d['seq'][:,0]==_m)[0]
                rdxdr[_m, :, :, :] += np.einsum('ijkl->jkl', d['rdxdr'][ids, :, :, :])
            d['rdxdr'] = rdxdr.reshape([N, L, 9])[:, :, [0, 4, 8, 1, 2, 5]]
            #d['rdxdr'] = np.einsum('ijklm->iklm', d['rdxdr'])\
            #.reshape([shp[0], shp[2], shp[3]*shp[4]])[:, :, [0, 4, 8, 1, 2, 5]]  #need to change
        
        #print(len(data['structure']))
        if self.base_potential:
            base_d = self.base_potential.calculate(data['structure'])
        else:
            base_d = {'energy': 0., 'force': 0., 'stress': 0.}
        
        #print(data['energy'])
        #print(base_d['energy'])

        d['energy'] = np.asarray(data['energy'] - base_d['energy'])
        d['force'] = np.asarray(data['force']) - base_d['force']
        if data['stress'] is not None:
            d['stress'] = np.asarray(data['stress']) - base_d['stress'] / units.GPa
        else:
            d['stress'] = data['stress'] 
        d['group'] = data['group']

        return d


def compute_descriptor(function, structure):
    """ Compute descriptor for one structure. """

    if function['type'] in ['BehlerParrinello', 'ACSF']:
        raise NotImplementedError

    elif function['type'] in ['wACSF', 'wacsf']:
        raise NotImplementedError

    elif function['type'] in ['SO4', 'Bispectrum', 'bispectrum']:
        raise NotImplementedError

    elif function['type'] in ['SO3', 'SOAP', 'soap']:
        d = DescriptorSO3(
                function['parameters']['nmax'],
                function['parameters']['lmax'],
                function['Rc'],
                derivative=True,
                stress=True
            ).calculate(structure)

    elif function['type'] in ['EAD', 'ead']:
        raise NotImplementedError

    elif function['type'] in ['SNAP', 'snap']:
        raise NotImplementedError
    else:
        msg = f"{function['type']} is not implemented"
        raise NotImplementedError(msg)
    
    # XXX ffr: What is this?
    if d['rdxdr'] is not None:
        N = d['x'].shape[0]
        L = d['x'].shape[1]
        rdxdr = np.zeros([N, L, 3, 3])
        for _m in range(N):
            ids = np.where(d['seq'][:,0]==_m)[0]
            rdxdr[_m, :, :, :] += np.einsum('ijkl->jkl', d['rdxdr'][ids, :, :, :])
        d['rdxdr'] = rdxdr.reshape([N, L, 9])[:, :, [0, 4, 8, 1, 2, 5]]
 
    return d




import ase.io
from ase import Atoms, Atom

def sort_atoms(atoms):
    atsymbs = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.pbc
    forces = atoms.get_forces()
    energy = atoms.get_potential_energy()

    idx_sorted = np.argsort(atsymbs)

    atoms_copy = Atoms()
    for i in idx_sorted:
        atoms_copy.append(Atom(atsymbs[i]))

    atoms_copy.set_positions(positions[idx_sorted])
    atoms_copy.set_cell(cell)
    atoms_copy.set_pbc(pbc)

    from ase.calculators.singlepoint import SinglePointCalculator
    stress = None # not yet needed
    magmoms = None
    calc = SinglePointCalculator(
        atoms_copy,
        energy=energy, forces=forces[idx_sorted]
    )
    atoms_copy.calc = calc
    
    return atoms_copy




# Modified by ffr, using ase.io, stress/virial is not yet read
# XXX For some part of the code it seems that we need to sort the atoms
def parse_xyz(structure_file):
    data = []
    for atoms in ase.io.iread(structure_file):
        atoms_p = sort_atoms(atoms)
        data.append({
            'structure': atoms_p,
            'energy': atoms_p.get_potential_energy(),
            'force': atoms_p.get_forces(),
            'stress': None,
            'group': 'random'
        })
    return data

def convert_to_ase_db(data, db_path='test.db'):
    from ase.db import connect
    
    with connect(db_path) as db:
        for d in data:
            struc = d['structure']
            d.pop('structure', None)
            db.write(struc, data=d)

def parse_ase_db(db_path, N=None, Random=False):
    from ase.db import connect

    data = []
    with connect(db_path) as db:
        if N is None:
            N = len(db)

        for i, row in enumerate(db.select()):
            structure = db.get_atoms(row.id)
            if "dft_stress" in row.data:
                stress = row.data["dft_stress"]
            elif "stress" in row.data:
                stress = row.data["stress"]
            else:
                stress = None

            if "group" in row.data:
                group = row.data["group"]
            else:
                group = None

            if "dft_energy" in row.data:
                eng = row.data["dft_energy"]
            else:
                eng = row.data["energy"]

            if "dft_force" in row.data:
                force = row.data["dft_force"]
            else:
                force = row.data["force"]
            data.append({'structure': structure,
                         'energy': eng,
                         'force': force,
                         'group': group,
                         'stress': stress})

            if i == (N-1):
                break

    return data

