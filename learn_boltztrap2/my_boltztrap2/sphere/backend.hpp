//    my_boltztrap2, a program for interpolating band structures and
//    calculating semi-classical transport coefficients.
//    Copyright (C) 2017-2025 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
//    Copyright (C) 2017-2025 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
//    Copyright (C) 2017-2025 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
//    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
//    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
//    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
//    Copyright (C) 2022 Robert Stanton <stantor@clarkson.edu>
//    Copyright (C) 2024 Haoyu (Daniel) Yang <yanghaoyu97@outlook.com>
//
//    This file is part of my_boltztrap2.
//
//    my_boltztrap2 is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    my_boltztrap2 is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with my_boltztrap2.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

#include <Eigen/Dense>

extern "C" {
#include <spglib.h>
}

/// Basic exception class used for all errors in this module
class Sphere_exception : public std::runtime_error
{
public:
    /// Basic constructor
    Sphere_exception(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

/// Type representing a lattice point with its squared norm attached
typedef std::pair<Eigen::Vector3i, double> point_with_sqnorm;


/// Different ways of expressing magnetic moments
/// unpolarized: no magnetization
/// collinear: collinear magnetization, one real number per atom
/// noncollinear: non-collinear magnetization, one 3-vector per atom
enum Magmom_type
{
    unpolarized = 0,
    collinear = 1,
    noncollinear = 2
};


/// Base class for classes that deal with symmetries in this module.
class Symmetry_analyzer
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] magmom - magnetic moments of each atom. It can be
    ///                     NULL for an unpolarized calculation, an
    ///                     array of natoms elements for collinear
    ///                     magnetism, or an array of 3 * natoms
    ///                     elements for noncollinear magnetism array of
    ///                     natoms elements. The kind of content is
    ///                     determined by mtype.
    /// @param[in] mtype - type of information in magmom (see
    ///                    description of the Magmom_type enum).
    /// @param[in] symprec - tolerance for symmetry search
    Symmetry_analyzer(double lattvec[3][3], double positions[][3],
                      const int types[], int natoms, double* magmom,
                      Magmom_type mtype, double symprec);
    /// @return the number of unique rotations.
    int get_nrotations() const
    {
        return static_cast<int>(this->rotations.size());
    }
    /// @return the cardinal of the basis of 3x3 tensors compatible with
    /// the symmetries.
    int get_ntensors()
    {
        if (!this->tensors_calculated) {
            this->compute_tensors();
        }
        return static_cast<int>(this->tensors.size());
    }
    /// Load the components of a 3x3 basis tensor into the argument.
    ///
    /// Here a basis is understood as a set of independent 3x3 tensor
    /// that span all possible 3x3 tensors compatible with the
    /// symmetries of the system. An order-2 linear response tensor must
    /// be a linear combination of these.
    /// Beware that the index is not checked
    /// @param[in] index - index of the point in the grid
    /// @param[out] tensor - array to store the components in row-major
    ///                      format
    void get_tensor(int index, double point[3])
    {
        if (!this->tensors_calculated) {
            this->compute_tensors();
        }
        for (int i = 0; i < 9; ++i) {
            point[i] = this->tensors[index](i / 3, i % 3);
        }
    }


protected:
    /// Set of rotations from the space group
    std::vector<Eigen::Matrix3i> rotations;
    /// Rotations in Cartesian coordinates
    std::vector<Eigen::Matrix3d> crotations;
    /// Analyze the structure with spglib and store the rotations of the
    /// space group, taking the magnetic configuration into account.
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] magmom - magnetic moments of each atom. It can be
    ///                     NULL for an unpolarized calculation, an
    ///                     array of natoms elements for collinear
    ///                     magnetism, or an array of 3 * natoms
    ///                     elements for noncollinear magnetism array of
    ///                     natoms elements. The kind of content is
    ///                     determined by mtype.
    /// @param[in] mtype - type of information in magmom (see
    ///                    description of the Magmom_type enum).
    /// @param[in] symprec - tolerance for symmetry search
    void analyze_symmetries(double lattvec[3][3], double positions[][3],
                            const int types[], int natoms,
                            double* magmom, Magmom_type mtype,
                            double symprec);

private:
    /// Has the basis of linear response tensors already been
    /// calculated?
    bool tensors_calculated;
    /// Once tensors_calculated is true, this vector will contain a
    /// basis for all 3x3 response tensors that are compatible with
    /// the symmetries.
    std::vector<Eigen::Matrix3d> tensors;
    /// Obtain a basis of 3x3 linear response tensors compatible with
    /// the symmetries of the system, put them in this->tensors and set
    /// this->tensors_calculated to true.
    void compute_tensors();
};

/// Base class for classes that build equivalences between points
class Equivalence_builder : public Symmetry_analyzer
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] magmom - magnetic moments of each atom. It can be
    ///                     NULL for an unpolarized calculation, an
    ///                     array of natoms elements for collinear
    ///                     magnetism, or an array of 3 * natoms
    ///                     elements for noncollinear magnetism array of
    ///                     natoms elements. The kind of content is
    ///                     determined by mtype.
    /// @param[in] mtype - type of information in magmom (see
    ///                    description of the Magmom_type enum).
    /// @param[in] symprec - tolerance for symmetry search
    Equivalence_builder(double lattvec[3][3], double positions[][3],
                        const int types[], int natoms, double* magmom,
                        Magmom_type mtype, double symprec)
        : Symmetry_analyzer(lattvec, positions, types, natoms, magmom,
                            mtype, symprec)
    {
    }
    /// Return a copy of the point mapping
    std::vector<int> get_mapping() const
    {
        return this->mapping;
    }

protected:
    /// Index of the equivalence class of each point of the grid
    std::vector<int> mapping;
};

/// Class used to classify equivalent lattice points in the intersection
/// between a supercell and a sphere
class Sphere_equivalence_builder : public Equivalence_builder
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] magmom - magnetic moments of each atom. It can be
    ///                     NULL for an unpolarized calculation, an
    ///                     array of natoms elements for collinear
    ///                     magnetism, or an array of 3 * natoms
    ///                     elements for noncollinear magnetism array of
    ///                     natoms elements. The kind of content is
    ///                     determined by mtype.
    /// @param[in] mtype - type of information in magmom (see
    ///                    description of the Magmom_type enum).
    /// @param[in] r - radius of the sphere
    /// @param[in] bounds - maximum absolute value of each direct
    /// coordinate
    /// in the supercell
    /// @param[in] symprec - tolerance for symmetry search
    Sphere_equivalence_builder(double lattvec[3][3],
                               double positions[][3], const int types[],
                               int natoms, double* magmom,
                               Magmom_type mtype, double r,
                               const int bounds[3], double symprec);
    /// Load the coordinates of a grid point into the argument
    ///
    /// Beware that the index is not checked
    /// @param[in] index - index of the point in the grid
    /// @param[out] point - array to store the coordinates of the point
    void get_point(int index, int point[3]) const
    {
        Eigen::Map<Eigen::Vector3i> wrapper(point);
        wrapper = this->grid[index].first;
    }

private:
    /// List of lattice points with their squared norm attached
    std::vector<point_with_sqnorm> grid;
    /// Fill the "grid" vector with all points in the intersection of
    /// the supercell and the sphere.
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] r - radius of the sphere
    /// @param[in] bounds - maximum absolute value of each direct
    /// coordinate in the supercell
    void create_grid(double lattvec[3][3], double r,
                     const int bounds[3]);
};

/// Class used to find all points in the BZ related to a single point by
/// symmetry.
class Degeneracy_counter : public Equivalence_builder
{
public:
    /// Basic constructor
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    /// @param[in] positions - atomic positions in direct coordinates
    /// @param[in] types - list of atom types
    /// @param[in] natoms - number of atoms in the structure
    /// @param[in] magmom - magnetic moments of each atom. It can be
    ///                     NULL for an unpolarized calculation, an
    ///                     array of natoms elements for collinear
    ///                     magnetism, or an array of 3 * natoms
    ///                     elements for noncollinear magnetism array of
    ///                     natoms elements. The kind of content is
    ///                     determined by mtype.
    /// @param[in] mtype - type of information in magmom (see
    ///                    description of the Magmom_type enum).
    /// @param[in] symprec - tolerance for symmetry search
    Degeneracy_counter(double lattvec[3][3], double positions[][3],
                       const int types[], int natoms, double* magmom,
                       Magmom_type mtype, double symprec);
    /// Obtain the number of symmetric versions of a point in reciprocal
    /// space, neglecting simple translations by a reciprocal lattice
    /// vector
    ///
    /// @param[in] point - point to analyze
    /// @return the degeneracy of the point
    int count_degeneracy(double point[3]);
    /// Load the coordinates of a point into the argument
    ///
    /// Beware that the index is not checked
    /// @param[in] index - index of the point in the grid
    /// @param[out] point - array to store the coordinates of the point
    void get_point(int index, double point[3]) const
    {
        for (int i = 0; i < 3; ++i) {
            point[i] = this->kpoints[index][i];
        }
    }

private:
    /// Tolerance for symmetry search
    double tolerance;
    /// Set of rotations from the space group, for vectors in reciprocal
    /// coordinates
    std::vector<Eigen::Matrix3d> krotations;
    /// List of reciprocal-space points
    std::vector<Eigen::Vector3d> kpoints;
};
