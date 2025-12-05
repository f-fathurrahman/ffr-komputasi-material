# -*- coding: utf-8 -*-
#    my_boltztrap2, a program for interpolating band structures and calculating
#    semi-classical transport coefficients.
#    Copyright (C) 2017-2025 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2025 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2025 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
#    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
#    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
#    Copyright (C) 2022 Robert Stanton <stantor@clarkson.edu>
#    Copyright (C) 2024 Haoyu (Daniel) Yang <yanghaoyu97@outlook.com>
#
#    This file is part of my_boltztrap2.
#
#    my_boltztrap2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    my_boltztrap2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with my_boltztrap2. If not, see <http://www.gnu.org/licenses/>.

# Implement support for loading and saving my_boltztrap2 calculations
# Calculations are saved in .bt2 files, which are simply xzip-compressed
# json files.

import copy
import datetime
import json
import lzma

import ase
import numpy as np

try:
    # Only needed for ASE >= 3.18.0
    import ase.cell

    _cell_object = True
except ModuleNotFoundError:
    _cell_object = False

import my_boltztrap2.dft
from my_boltztrap2.misc import warning
from my_boltztrap2.version import *


def atoms2dict(atoms):
    """Convert an ase Atoms object into a dictionary of json-izable objects.

    Note that this is not a general-purpose implementation. Featured not used
    by my_boltztrap2 have not been tested.
    """
    nruter = dict()
    nruter["BoltzTraP2_type"] = "Atoms"
    nruter["cell"] = atoms.cell
    nruter["pbc"] = atoms.pbc
    nruter["atoms"] = []
    for a in atoms:
        try:
            magmom = a.magmom
        except AttributeError:
            magmom = a.magmom
        atom = dict(
            symbol=a.symbol,
            position=a.position,
            tag=a.tag,
            momentum=a.momentum,
            mass=a.mass,
            magmom=magmom,
            charge=a.charge,
        )
        nruter["atoms"].append(atom)
    return nruter


def dict2atoms(d):
    """Convert a dictionary loaded from a json file into an Atoms object."""
    atoms = []
    for a in d["atoms"]:
        atoms.append(ase.Atom(**a))
    return ase.Atoms(atoms, cell=d["cell"], pbc=d["pbc"])


def data2dict(data):
    """Convert a DFTData object into a a dictionary of json-izable objects."""
    nruter = dict()
    nruter["BoltzTraP2_type"] = "DFTData"
    nruter["ebands"] = data.ebands
    nruter["mommat"] = data.mommat
    nruter["kpoints"] = data.kpoints
    nruter["atoms"] = data.atoms
    nruter["fermi"] = data.fermi
    nruter["dosweight"] = data.dosweight
    nruter["nelect"] = data.nelect
    nruter["magmom"] = data.magmom
    return nruter


class BT2Loader:
    """Loader of bt2 files."""

    def __init__(self, d):
        # This is a loader in the style of those in the DFT module, but
        # taking a dictionary as its argument.
        if not isinstance(d, dict):
            raise my_boltztrap2.dft.LoaderError(
                "this loader only works with dictionaries"
            )
        try:
            self.ebands = np.array(d["ebands"])
            self.mommat = d["mommat"]
            if self.mommat is not None:
                self.mommat = np.array(self.mommat)
            self.kpoints = np.array(d["kpoints"])
            self.atoms = d["atoms"]
            self.sysname = self.atoms.get_chemical_formula()
            self.fermi = d["fermi"]
            self.dosweight = d["dosweight"]
            self.nelect = d["nelect"]
            self.magmom = d["magmom"]
        except KeyError:
            raise my_boltztrap2.dft.LoaderError("a key is missing from the input")


my_boltztrap2.dft.register_loader("BT2", BT2Loader)


def dict2data(d):
    """Convert a dictionary loaded from a json file into a DFTData object."""
    return my_boltztrap2.dft.DFTData(d, derivatives=True)


def object_hook(d):
    """Build a my_boltztrap2-related object from d when possible.

    If d is a dictionary that looks like it can be converted to some kind of
    object supported by this module, perform the conversion and return the
    result. Otherwise, just return d unchanged.
    """
    if "BoltzTraP2_type" in d:
        if d["BoltzTraP2_type"] == "Array":
            if d["array_type"] == "Complex":
                return np.array(d["real"]) + 1j * np.array(d["imag"])
            else:
                return np.array(d["data"])
        if d["BoltzTraP2_type"] == "Atoms":
            return dict2atoms(d)
        if d["BoltzTraP2_type"] == "DFTData":
            return dict2data(d)
        return d
    return d


class JSONEncoder(json.JSONEncoder):
    """Class that extends JSONEncoder to handle BoltzTrap2-related objects."""

    def default(self, o):
        """Return a json-izable version of o or delegate on the base class."""
        if isinstance(o, np.generic):
            # Deal with non-serializable types such as numpy.int64
            return o.item()
        if isinstance(o, np.ndarray):
            # Arrays of complex numbers are serialized as a dictionary.
            nruter = dict(BoltzTraP2_type="Array")
            if np.iscomplexobj(o):
                nruter["array_type"] = "Complex"
                nruter["real"] = o.real.tolist()
                nruter["imag"] = o.imag.tolist()
            # Other arrays are simply serialized as a list with a type
            else:
                nruter["array_type"] = "NumPy/" + o.dtype.name
                nruter["data"] = o.tolist()
            return nruter
        if _cell_object:
            if isinstance(o, ase.cell.Cell):
                return o.array
        if isinstance(o, ase.Atoms):
            return atoms2dict(o)
        if isinstance(o, my_boltztrap2.dft.DFTData):
            return data2dict(o)
        return json.JSONEncoder.default(self, o)


def save_calculation(filename, data, equivalences, coeffs, metadata):
    """Save the objects describing a calculation to a xzipped json file.

    Args:
        filename: path to the file
        data: DFTData object
        equivalences: list of k-point equivalence classes
        coeffs: interpolation coefficients
        metadata: metadata dictionary

    Returns:
        None
    """
    with lzma.open(filename, "wt", encoding="utf-8", preset=0) as f:
        # Calling dumps() instead of dump() makes this function twice as
        # fast at the cost of requiring more memory. Currently, the
        # implementation of dumps() is pure Python whereas .dump() is
        # optimized in C.
        f.write(
            json.dumps([data, equivalences, coeffs, metadata], cls=JSONEncoder)
        )


def load_calculation(filename):
    """Load the objects describing a calculation from a xzipped json file.

    Args:
        filename: path to the file to be read

    Returns:
        A sequence of four elements:
        1. A DFTData object
        2. A list of k-point equivalence classes
        3. The interpolation coefficients
        4. A metadata dictionary
    """
    with lzma.open(filename, "rt", encoding="utf-8") as f:
        nruter = json.load(f, object_hook=object_hook)
    nruter[1] = [np.array(i) for i in nruter[1]]
    nruter[2] = np.array(nruter[2])
    if nruter[3]["format_version"] != BT2_FORMAT_VERSION:
        warning(
            "unexpected bt2 format version ({} instead of {})".format(
                nruter[3]["format_version"], BT2_FORMAT_VERSION
            )
        )
    return nruter


def save_results(
    filename,
    data,
    fermi,
    Tr,
    mu0,
    mur,
    N,
    sdos,
    cv,
    cond,
    seebeck,
    kappa,
    hall,
    metadata,
):
    """Save all the results from an integration step to a single compressed
    JSON file.

    Args:
        filename: path to the file
        data: DFTData object
        fermi: Fermi level
        Tr: vector of temperatures
        mu0: vector of intrinsic chemical potentials at each temperatures
        mur: vector of chemical potentials
        N: results for electron count
        sdos: results for the smoothed DOS
        cv: results for the heat capacity
        cond: results for the electric conductivity
        seebeck: results for the Seebeck coefficients
        kappa: results for the carrier contribution to the thermal conductivity
        hall: results for the Hall coefficient
        metadata: metadata dictionary

    Returns:
        None.
    """
    with lzma.open(filename, "wt", encoding="utf-8", preset=0) as f:
        f.write(
            json.dumps(
                [
                    data,
                    fermi,
                    Tr,
                    mu0,
                    mur,
                    N,
                    sdos,
                    cv,
                    cond,
                    seebeck,
                    kappa,
                    hall,
                    metadata,
                ],
                cls=JSONEncoder,
            )
        )


def load_results(filename):
    """Load the results from an integration step stored in a compressed
    JSON file.

    Args:
        filename: path to the file to be loaded

    Returns:
        A sequence of 13 elements:
        1. DFTData object
        2. Fermi level
        3. vector of temperatures
        4. vector of intrinsic chemical potentials at each temperatures
        5. vector of chemical potentials
        6. results for electron count
        7. results for the smoothed DOS
        8. results for the heat capacity
        9. results for the electric conductivity
        10. results for the Seebeck coefficients
        11. results for the carrier contribution to the thermal conductivity
        12. results for the Hall coefficient
        13. metadata dictionary
    """
    with lzma.open(filename, "rt", encoding="utf-8") as f:
        nruter = json.load(f, object_hook=object_hook)
    for i in range(2, len(nruter) - 1):
        nruter[i] = np.array(nruter[i])
    if nruter[12]["format_version"] != BTJ_FORMAT_VERSION:
        warning(
            "unexpected btj format version ({} instead of {})".format(
                nruter[12]["format_version"], BTJ_FORMAT_VERSION
            )
        )
    return nruter


def load_basic(filename):
    """Load the 'data' and 'metadata' objects from a .bt2 or a .btj file.

    Args:
        filename: path to the file to be loaded

    Returns:
        A sequence of two elements: a DFTData object and a metadata dictionary.
    """
    with lzma.open(filename, "rt", encoding="utf-8") as f:
        nruter = json.load(f, object_hook=object_hook)
    return (nruter[0], nruter[-1])


def gen_bt2_metadata(data, derivatives):
    """Create a metadata dictionary for a bt2 file.

    Args:
        data: DFTData object describing the calculation
        derivatives: logical variable stating whether the calculation included
            derivative information in the fit.

    Returns:
        A dictionary of metadata.
    """
    return dict(
        program_version=PROGRAM_VERSION,
        format_version=BT2_FORMAT_VERSION,
        ctime=datetime.datetime.now().isoformat(),
        source=data.source,
        derivatives=derivatives,
    )


def gen_btj_metadata(metadata0, scattering_model):
    """Create a metadata dictionary for a btj file.

    Args:
        metadata0: metadata dictionary from the input bt2 file
        scattering_model: string or array with the scattering model used in the
            calculation

    Returns:
        A new dictionary of metadata.
    """
    metadata = copy.copy(metadata0)
    metadata["format_version"] = BTJ_FORMAT_VERSION
    metadata["ctime"] = datetime.datetime.now().isoformat()
    if isinstance(scattering_model, str):
        sm = scattering_model
    else:
        sm = "custom"
    metadata["scattering_model"] = sm
    return metadata
