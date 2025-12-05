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


#include <cmath>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <set>
#include <utility>

#include "backend.hpp"

/// Simple template for a sign function
/// @param[in] x - an argument for which "<" works and that can be
///                constructed passing 0 as the argument
/// @return -1, 0 or 1 depending on the sign of the argument
template<typename T> int sign(T x)
{
    return (T(0) < x) - (x < T(0));
}

/// Find a permutation between two sets of atomic positions.
///
/// The coordinates are assumed to be direct, and translations by an
/// integer are ignored. If the permutation does not exist, throw a
/// Sphere_exception.
/// @param[in] reference - 3 x natoms array with the original direct
///                        coordinates.
/// @param[in] transformed - 3 x natoms array with the transformed direct
///                          coordinates.
/// @param[in] symprec - tolerance for symmetry search.
/// @return - the permutation as an array of indices.
std::vector<std::size_t> find_permutation(
    const Eigen::Ref<const Eigen::MatrixXd>& reference,
    const Eigen::Ref<const Eigen::MatrixXd>& transformed,
    double symprec)
{
    if (reference.rows() != 3 || transformed.rows() != 3) {
        throw Sphere_exception("inputs must be n x 3 arrays");
    }
    auto natoms = reference.cols();
    if (transformed.cols() != natoms) {
        throw Sphere_exception("reference and transformed positions"
                               " must have the same shape");
    }
    std::vector<std::size_t> nruter;
    std::vector<bool> available(natoms);
    std::fill(available.begin(), available.end(), true);
    for (int ia = 0; ia < natoms; ++ia) {
        int minr;
        int minc;
        Eigen::ArrayXd delta(natoms);
        for (int ja = 0; ja < natoms; ++ja) {
            if (available[ja]) {
                Eigen::Array3d diff =
                    (reference.col(ja) - transformed.col(ia)).array();
                diff = (diff - diff.round()).abs();
                delta(ja) = diff.sum();
            }
            else {
                delta(ja) = std::numeric_limits<double>::infinity();
            }
        }
        double min = delta.minCoeff(&minr, &minc);
        if (min >= symprec) {
            throw Sphere_exception("cannot translate a symmetry into"
                                   " an atom permutation");
        }
        available[minr] = false;
        nruter.emplace_back(minr);
    }
    return nruter;
}

/// Find the permutation associated to each symmetry operation.
///
/// @param[in] lattvec - lattice vectors, as columns
/// @param[in] positions - atomic positions in direct coordinates
/// @param[in] types - list of atom types
/// @param[in] natoms - number of atoms in the structure
/// @param[in] data - information about the symmetry operations of the system
/// @return a vector of permutations, one for each operation.
std::vector<std::vector<std::size_t>> get_permutations(
    double lattvec[3][3], double positions[][3], const int types[],
    int natoms, SpglibDataset* data, double symprec)
{
    int nops = data->n_operations;
    std::vector<std::vector<std::size_t>> nruter;
    nruter.reserve(nops);
    Eigen::Map<Eigen::MatrixXd> direct(&(positions[0][0]), 3, natoms);
    for (int iop = 0; iop < nops; ++iop) {
        Eigen::Map<Eigen::Matrix3i> rwrapper(
            &(data->rotations[iop][0][0]));
        Eigen::Map<Eigen::Vector3d> twrapper(
            &(data->translations[iop][0]));
        Eigen::MatrixXd transformed =
            rwrapper.cast<double>().transpose() * direct;
        for (int iatom = 0; iatom < natoms; ++iatom) {
            transformed.col(iatom) += twrapper;
        }
        nruter.emplace_back(
            find_permutation(direct, transformed, symprec));
    }
    return nruter;
}

/// Objects of this simple POD class describe the effect of a symmetry
/// operation on the atomic magnetic moments
class Compatibility_verdict
{
public:
    /// Default constructorl
    Compatibility_verdict() : forward(false), backward(false)
    {
    }
    /// Can the operation preserve the atomic magnetizations without
    /// reversing time?
    bool forward;
    /// Can the operation preserve the atomic magnetizations after
    /// reversing time?
    bool backward;
};

/// Determine whether a symmetry operation is compatible with the
/// magnetic configuration of the system.
///
/// @param[in] rotation - rotation matrix associated to the symmetry
/// @param[in] permutation - atomic permutation associated to the
///                          symmetry
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
/// @return true or false, depending on whether
Compatibility_verdict determine_compatibility(
    const Eigen::Ref<const Eigen::Matrix3d>& rotation,
    const std::vector<std::size_t>& permutation, double* magmom,
    Magmom_type mtype, double symprec)
{
    Compatibility_verdict nruter;
    if (mtype == Magmom_type::unpolarized) {
        // Short circuit for unpolarized calculations
        nruter.forward = true;
        nruter.backward = true;
    }
    else if (mtype == Magmom_type::collinear) {
        // For collinear calculations, a global inversion of the spins
        // can be performed independently of the lattice. Hence time
        // reversal needs not be considered.
        auto natoms = static_cast<int>(permutation.size());
        Eigen::Map<Eigen::ArrayXd> oldm(magmom, natoms);
        Eigen::ArrayXd newm(natoms);
        for (int ia = 0; ia < natoms; ++ia) {
            newm(ia) = oldm(permutation[ia]);
        }
        nruter.forward =
            std::min((oldm - newm).abs().maxCoeff(),
                     (oldm + newm).abs().maxCoeff()) < symprec;
        nruter.backward = nruter.forward;
    }
    else {
        // Case with noncollinear magnetization. The magnetic moments
        // are rotated, and time reversal needs to be taken into
        //  account.
        auto natoms = static_cast<int>(permutation.size());
        // magmom will be a flat representation of a row-major array of
        // size natoms x 3, which we interpret as a column-major array
        // of size 3 x natoms.
        Eigen::Map<
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
            oldm(magmom, 3, natoms);
        Eigen::MatrixXd newm(3, natoms);
        for (int ia = 0; ia < natoms; ++ia) {
            newm.col(ia) = oldm.col(permutation[ia]);
        }
        // Transform each magnetic moment according to the symmetries
        // If the operation is an improper rotation, there is an extra
        // change of sign involved, since magnetic moments are axial
        // vectors.
        double prefactor =
            static_cast<float>(sign(rotation.determinant()));
        newm = (prefactor * rotation * newm).eval();

        if ((oldm - newm).norm() < symprec) {
            nruter.forward = true;
            nruter.backward = false;
        }
        else if ((oldm + newm).norm() < symprec) {
            nruter.forward = false;
            // TODO: this should probably be true, but does VASP agree?
            nruter.backward = false;
        }
    }
    return nruter;
}

Symmetry_analyzer::Symmetry_analyzer(double lattvec[3][3],
                                     double positions[][3],
                                     const int types[], int natoms,
                                     double* magmom, Magmom_type mtype,
                                     double symprec)
    : tensors_calculated(false)
{
    // Obtain the rotation matrices, taking the magnetic configuration
    // into account.
    this->analyze_symmetries(lattvec, positions, types, natoms, magmom,
                             mtype, symprec);
}

/// Trivial "less than" operator implementation for Eigen::Matrix3i
namespace std
{
template<> struct less<Eigen::Matrix3i>
{
    bool operator()(const Eigen::Ref<const Eigen::Matrix3i>& a,
                    const Eigen::Ref<const Eigen::Matrix3i>& b) const
    {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (a(i, j) < b(i, j)) {
                    return true;
                }
                if (a(i, j) > b(i, j)) {
                    return false;
                }
            }
        }
        return false;
    }
};
} // namespace std

void Symmetry_analyzer::analyze_symmetries(
    double lattvec[3][3], double positions[][3], const int types[],
    int natoms, double* magmom, Magmom_type mtype, double symprec)
{
    // Call spglib to get all information about symmetries
    SpglibDataset* data =
        spg_get_dataset(lattvec, positions, types, natoms, symprec);
    if (data == NULL) {
        throw Sphere_exception("spglib's spg_get_dataset"
                               " returned a NULL pointer");
    }
    // Obtain the atomic permutations associated to each symmetry
    // operation.
    auto permutations = get_permutations(lattvec, positions, types,
                                         natoms, data, symprec);
    // Prepare a quick way to change between Cartesian and direct
    // coordinates.
    Eigen::Map<Eigen::Matrix3d> tlattvec(&(lattvec[0][0]));
    Eigen::PartialPivLU<Eigen::Matrix3d> solver(tlattvec);
    // Scan through the operations and keep only a unique set of
    // rotations plus the corresponding roto-inversions. Discard those
    // that are incompatible with the magnetic configuration.
    int nops = data->n_operations;
    std::set<Eigen::Matrix3i> unique;
    for (int i = 0; i < nops; ++i) {
        Eigen::Map<Eigen::Matrix3i> wrapper(
            &(data->rotations[i][0][0]));
        Eigen::Matrix3d crotation =
            solver.solve(wrapper.cast<double>() * tlattvec).transpose();
        auto compatibility = determine_compatibility(
            crotation, permutations[i], magmom, mtype, symprec);
        if (compatibility.forward) {
            unique.insert(wrapper.transpose());
        }
        if (compatibility.backward) {
            unique.insert(-wrapper.transpose());
        }
    }
    this->rotations.assign(unique.begin(), unique.end());
    nops = static_cast<int>(this->rotations.size());
    // Fill the crotations vector with the Cartesian rotations
    for (int i = 0; i < nops; ++i) {
        Eigen::Matrix3d crotation =
            solver
                .solve(this->rotations[i].cast<double>().transpose() *
                       tlattvec)
                .transpose();
        this->crotations.emplace_back(crotation);
    }
    // Free up the space allocated by spglib
    spg_free_dataset(data);
}

void Symmetry_analyzer::compute_tensors()
{
    // Linear response tensors must be invariant with respect to a
    // conjugation by any of the rotation matrices.
    // Express these constraints as the coefficients of a linear system
    // where the unknown (a generic 3x3 tensor) is raveled to a 9-vector
    // in C order.
    Eigen::Matrix3d one(Eigen::Matrix3d::Identity());
    auto nops = this->rotations.size();
    Eigen::MatrixXd constraints(9 * nops + 3, 9);
    for (std::size_t iop = 0; iop < nops; ++iop) {
        // The constraint matrix is based on the Cartesian
        // representation of the rotation.
        Eigen::Matrix3d rot = this->crotations[iop];
        auto offset = 9 * iop;
        for (int alpha = 0; alpha < 3; ++alpha) {
            for (int delta = 0; delta < 3; ++delta) {
                auto row = offset + 3 * alpha + delta;
                for (int beta = 0; beta < 3; ++beta) {
                    for (int gamma = 0; gamma < 3; ++gamma) {
                        auto col = 3 * beta + gamma;
                        constraints(row, col) =
                            rot(beta, alpha) * rot(gamma, delta) -
                            one(beta, alpha) * one(gamma, delta);
                    }
                }
            }
        }
    }
    // Insert additional constraints to ensure that the solution are
    // symmetric tensors regardless of the space group.
    constraints.bottomRows<3>().fill(0.);
    constraints(9 * nops, 3 * 0 + 1) = 1.;
    constraints(9 * nops, 3 * 1 + 0) = -1.;
    constraints(9 * nops + 1, 3 * 0 + 2) = 1.;
    constraints(9 * nops + 1, 3 * 2 + 0) = -1.;
    constraints(9 * nops + 2, 3 * 1 + 2) = 1.;
    constraints(9 * nops + 2, 3 * 2 + 1) = -1.;

    // Use an LU decomposition to extract a basis for the kernel of the
    // coefficient matrix, and recast each element of the basis as a
    // 3x3 matrix.
    Eigen::FullPivLU<Eigen::MatrixXd> lu(constraints);
    Eigen::MatrixXd kernel = lu.kernel();
    auto nbasis = kernel.cols();
    for (int i = 0; i < nbasis; ++i) {
        Eigen::Matrix3d tensor;
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                tensor(r, c) = kernel(3 * r + c, i);
            }
        }
        this->tensors.emplace_back(tensor);
    }

    this->tensors_calculated = true;
}

/// Simple class for computing squared norms of lattice vectors
class Lattice_tape
{
public:
    /// Basic constructor.
    ///
    /// @param[in] lattvec - lattice vectors, as columns
    Lattice_tape(const double lattvec[3][3])
    {
        Eigen::Map<const Eigen::Matrix3d> wrapper(&(lattvec[0][0]));
        this->metric = wrapper * wrapper.transpose();
    }
    /// Compute the squared norm of a lattice vector
    ///
    /// @param[in] v - vector in direct coordinates
    /// @return the squared norm of v
    double measure(const Eigen::Ref<const Eigen::Vector3i>& v) const
    {
        Eigen::Vector3d dv = v.cast<double>();
        return dv.transpose() * this->metric * dv;
    }

private:
    /// Metric tensor
    Eigen::Matrix3d metric;
};

/// Comparison function for pairs or tuples based on the second element
///
/// @param p1 - first operand
/// @param p2 - second operand
/// @return the result of p1.second < p2.second
template<typename T> bool compare_second(const T& p1, const T& p2)
{
    return p1.second < p2.second;
}

Sphere_equivalence_builder::Sphere_equivalence_builder(
    double lattvec[3][3], double positions[][3], const int types[],
    int natoms, double* magmom, Magmom_type mtype, double r,
    const int bounds[3], double symprec)
    : Equivalence_builder(lattvec, positions, types, natoms, magmom,
                          mtype, symprec)
{
    // Generate a list of all lattice points in the intersection and
    // store their squared norms for later use
    this->create_grid(lattvec, r, bounds);
    // Initialize the point mapping
    this->mapping.resize(this->grid.size());
    std::fill(this->mapping.begin(), this->mapping.end(), -1);
    // Find an equivalence class for each point
    point_with_sqnorm lo = std::make_pair(Eigen::Vector3i::Zero(), 0.);
    point_with_sqnorm hi = std::make_pair(Eigen::Vector3i::Zero(), 0.);
    for (std::size_t i = 0; i < this->grid.size(); ++i) {
        // Each point can only be equivalent to a point with
        // the same norm. Some room is left for rounding errors
        Eigen::Vector3i point = this->grid[i].first;
        double sqnorm = this->grid[i].second;
        lo.second = (1. - symprec) * sqnorm;
        hi.second = (1. + symprec) * sqnorm;
        std::size_t lbound = std::distance(
            this->grid.begin(),
            std::lower_bound(this->grid.begin(), this->grid.end(), lo,
                             compare_second<point_with_sqnorm>));
        std::size_t ubound = std::distance(
            this->grid.begin(),
            std::upper_bound(this->grid.begin(), this->grid.end(), hi,
                             compare_second<point_with_sqnorm>));
        ubound = std::min(ubound, i);
        for (std::size_t o = 0; o < this->rotations.size(); ++o) {
            // Apply each symmetry operation to the point
            Eigen::Vector3i image = this->rotations[o] * point;
            // And compare the result to each of the candidates
            for (std::size_t j = lbound; j < ubound; ++j) {
                if (image == this->grid[j].first) {
                    // If they match, find the representative of the
                    // equivalence class, assign this point to the class
                    // and exit the loop
                    auto k = j;
                    while (this->mapping[k] != k) {
                        k = this->mapping[k];
                    }
                    this->mapping[i] = static_cast<int>(k);
                    break;
                }
            }
            if (this->mapping[i] != -1) {
                break;
            }
        }
        // If no equivalence class was found for this point, start a
        // new one.
        if (this->mapping[i] == -1) {
            this->mapping[i] = static_cast<int>(i);
        }
    }
}

void Sphere_equivalence_builder::create_grid(double lattvec[3][3],
                                             double r,
                                             const int bounds[3])
{
    double r2(r * r);
    Lattice_tape tape(lattvec);
    for (int i = -bounds[0]; i <= bounds[0]; ++i) {
        for (int j = -bounds[1]; j <= bounds[1]; ++j) {
            for (int k = -bounds[2]; k <= bounds[2]; ++k) {
                Eigen::Vector3i point(i, j, k);
                double n2 = tape.measure(point);
                if (n2 < r2) {
                    this->grid.push_back(std::make_pair(point, n2));
                }
            }
        }
    }
    // The list is stored in order of increasing norm
    std::sort(this->grid.begin(), this->grid.end(),
              compare_second<point_with_sqnorm>);
}

Degeneracy_counter::Degeneracy_counter(
    double lattvec[3][3], double positions[][3], const int types[],
    int natoms, double* magmom, Magmom_type mtype, double symprec)
    : Equivalence_builder(lattvec, positions, types, natoms, magmom,
                          mtype, symprec),
      tolerance(symprec)
{
    // Obtain the rotation matrices in the reciprocal basis
    Eigen::Map<Eigen::Matrix3d> tlattvec(&(lattvec[0][0]));
    Eigen::Matrix3d metric = tlattvec * tlattvec.transpose();
    Eigen::ColPivHouseholderQR<Eigen::Matrix3d> solver;
    solver.compute(metric);
    for (std::size_t ir = 0; ir < this->rotations.size(); ++ir) {
        this->krotations.push_back(
            (solver.solve(
                 this->rotations[ir].cast<double>().transpose() *
                 metric))
                .transpose());
    }
}

int Degeneracy_counter::count_degeneracy(double point[3])
{
    // Fill the list with all possible images of the point
    this->kpoints.clear();
    Eigen::Map<Eigen::Vector3d> kpoint(point);
    for (std::size_t o = 0; o < this->krotations.size(); ++o) {
        Eigen::Vector3d image = this->krotations[o] * kpoint;
        image -= image.array().round().matrix().eval();
        this->kpoints.push_back(image);
    }
    // Initialize the point mapping
    this->mapping.resize(this->kpoints.size());
    std::fill(this->mapping.begin(), this->mapping.end(), -1);
    // Build the mapping by finding points related by translations
    for (std::size_t i = 0; i < this->kpoints.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            Eigen::Vector3d delta = this->kpoints[i] - this->kpoints[j];
            delta -= delta.array().round().matrix().eval();
            if (delta.cwiseAbs().maxCoeff() < this->tolerance) {
                auto k = j;
                while (this->mapping[k] != k) {
                    k = this->mapping[k];
                }
                this->mapping[i] = static_cast<int>(k);
                break;
            }
        }
        if (this->mapping[i] == -1) {
            this->mapping[i] = static_cast<int>(i);
        }
    }
    /// Count and return the number of classes in the mapping
    std::sort(this->mapping.begin(), this->mapping.end());
    return static_cast<int>(std::unique(this->mapping.begin(), this->mapping.end()) -
           this->mapping.begin());
}
