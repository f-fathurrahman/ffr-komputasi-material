/* Copyright 2016 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** @package MyTBTKcalc
 *  @file PropertyExtractor.h
 *  @brief Base class for PropertyExtractors.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_PROPERTY_EXTRACTOR
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_PROPERTY_EXTRACTOR

#include "MyTBTK/HoppingAmplitudeSet.h"
#include "MyTBTK/Index.h"
#include "MyTBTK/Property/AbstractProperty.h"
#include "MyTBTK/Property/Density.h"
#include "MyTBTK/Property/DOS.h"
#include "MyTBTK/Property/LDOS.h"
#include "MyTBTK/Property/Magnetization.h"
#include "MyTBTK/Property/SpinPolarizedLDOS.h"
#include "MyTBTK/Solver/Solver.h"

#include <complex>
//#include <initializer_list>

namespace MyTBTK{
namespace PropertyExtractor{

/** @brief Base class for PropertyExtractors.
 *
 *  A PropertyExtractor can be used to extract @link Property::AbstractProperty
 *  Properties@endlink from a @link Solver::Solver Solver@endlink. The
 *  PropertyExtractor::PropertyExtractor is a base class for such
 *  @link PropertyExtractor PropertyExtractors@endlink, which each corresponds
 *  to a particular @link Solver::Solver Solver@endlink. See the documentation
 *  for Diagonalizer, BlockDiagonalizer, ArnoldiIterator, and ChebyshevExpander
 *  for examples of specific production ready @link PropertyExtractor
 *  PropertyExtractors@endlink.
 *
 *  The @link PropertyExtractor::PropertyExtractor PropertyExtractors@endlink
 *  provide a uniform interface to @link Solver::Solver Solvers@endlink and
 *  allow for @link Property::AbstractProperty Properties@endlink to be
 *  extracted with limited knowledge about @link Solver::Solver Solver@endlink
 *  specific details. The use of @link
 *  PropertyExtractor::AbstractPropertyExtractor PropertyExtractors@endlink
 *  also makes it possible to switch between different @link Solver::Solver
 *  Solvers@endlink with minimal changes to the code.
 *
 *  # Example
 *  \snippet PropertyExtractor/PropertyExtractor.cpp PropertyExtractor */
class PropertyExtractor{
public:
	/** Constructs a PropertyExtractor::PropertyExtractor. */
	PropertyExtractor();

	/** Destructor. */
	virtual ~PropertyExtractor();

	/** Set the energy window used for energy dependent quantities. The
	 *  energy window is set to be real.
	 *
	 *  @param lowerBound The lower bound for the energy window.
	 *  @param upperBound The upper bound for the energy window.
	 *  @param energyResolution The number of energy points used to resolve
	 *  the energy window. */
	virtual void setEnergyWindow(
		double lowerBound,
		double upperBound,
		int energyResolution
	);

	/** Set the energy window used for energy dependent quantities. The
	 *  energy window is set to consist of Matsubara energies.
	 *
	 *  @param lowerFermionicMatsubaraEnergyIndex The lower Fermionic
	 *  Matsubara energy index.
	 *
	 *  @param upperFermionicMatsubaraEnergyIndex The upper Fermionic
	 *  Matsubara energy index.
	 *
	 *  @param lowerBosonicMatsubaraEnergyIndex The lower Bosonic
	 *  Matsubara energy index.
	 *
	 *  @param upperBosonicMatsubaraEnergyIndex The upper Bosonic
	 *  Matsubara energy index. */
	virtual void setEnergyWindow(
		int lowerFermionicMatsubaraEnergyIndex,
		int upperFermionicMatsubaraEnergyIndex,
		int lowerBosonicMatsubaraEnergyIndex,
		int upperBosonicMatsubaraEnergyIndex
	);

	/** Set the size of the energy infinitesimal that can be used to add
	 *  for example an \f$i\delta\f$ term to the denominator of the Green's
	 *  function.
	 *
	 *  @param energyInfinitesimal The energy infinitesimal \f$\delta\f$.
	 */
	virtual void setEnergyInfinitesimal(double energyInfinitesimal);

	/** Calculate the density. This function should be overriden by those
	 *  deriving classes that provide support for calculating the density.
	 *  By default the PropertyExtractor prints an error message that the
	 *  given property is not supported.
	 *
	 *  @param pattern Specifies the index pattern for which to calculate
	 *  the density. For example, assume that the index scheme is
	 *  {x, y, z, spin}. {ID_X, 5, 10, IDX_SUM_ALL} will calculate the
	 *  density for each x along (y,z)=(5,10) by summing over spin.
	 *  Similarly {ID_X, 5, IDX_Y, IDX_SUM_ALL} will return a two
	 *  dimensional density for all x and z and y = 5. Note that IDX_X
	 *  IDX_Y, and IDX_Z refers to the first, second, and third index used
	 *  by the routine to create a one-, two-, or three-dimensional output,
	 *  rather than being tied to the x, y, and z used as physical
	 *  subindices.
	 *
	 *  @param ranges Speifies the number of elements for each subindex. Is
	 *   ignored for indices specified with positive integers in the
	 *  pattern, but is used to loop from 0 to the value in ranges for
	 *  IDX_X, IDX_Y, IDX_Z, and IDX_SUM_ALL. Appropriate ranges
	 *  corresponding to the two pattern examples above are
	 *  {SIZE_X, 1, 1, NUM_SPINS} and {SIZE_X, 1, SIZE_Z, NUM_SPINS},
	 *  respectively.
	 *
	 *  @return A density array with size equal to the number of points
	 *  included by specified patter-range combination. */
	virtual Property::Density calculateDensity(
		Index pattern,
		Index ranges
	);

	/** Calculate the density. This function should be overriden by those
	 *  deriving classes that provide support for calculating the density.
	 *  By default the PropertyExtractor prints an error message that the
	 *  given property is not supported.
	 *
	 *  @param patterns A list of patterns that will be matched against the
	 *  @link Index Indices @endlink in the Model to determine which @link
	 *  Index Indices @endlink for which to calculate the Density.
	 *
	 *  @return A Property::Density for the @link Index Indices @endlink
	 *  that match the patterns. */
	virtual Property::Density calculateDensity(
		std::vector<Index> patterns
	);

	/** Calculate the magnetization. This function should be overriden by
	 *  those deriving classes that provide support for calculating the
	 *  magnetization. By default the PropertyExtractor prints an error
	 *  message that the given property is not supported.
	 *
	 *  @param pattern Specifies the index pattern for which to calculate
	 *  the magnetization. For example, assume that the index scheme is
	 *  {x, y, z, spin}. {ID_X, 5, 10, IDX_SPIN} will calculate the
	 *  magnetization for each x along (y,z)=(5,10). Similarly
	 *  {ID_X, 5, IDX_Y, IDX_SPIN} will return a two dimensional
	 *  magnetiation for all x and z and y = 5. Note that IDX_X, IDX_Y, and
	 *  IDX_Z refers to the first, second, and third index used by the
	 *  routine to create a one-, two-, or three-dimensional output, rather
	 *  than being tied to the x, y, and z used as physical subindices.
	 *
	 *  @param ranges Speifies the number of elements for each subindex. Is
	 *  ignored for indices specified with positive integers in the
	 *  pattern, but is used to loop from 0 to the value in ranges for
	 *  IDX_X, IDX_Y, IDX_Z, and IDX_SUM_ALL. Appropriate ranges
	 *  corresponding to the two pattern examples above are
	 *  {SIZE_X, 1, 1, NUM_SPINS} and {SIZE_X, 1, SIZE_Z, NUM_SPINS},
	 *  respectively.
	 *
	 *  @return A magnetization array with size equal to four times the
	 *  number of points included by specified patter-range combination.
	 *  The four entries are
	 *  \f[
	 *      \left[\begin{array}{cc}
	 *          0   & 1\\
	 *          2   & 3
	 *      \end{array}\right] =
	 *      \left[\begin{array}{cc}
	 *          \langle c_{i\uparrow}^{\dagger}c_{i\uparrow}\rangle         & \langle c_{i\uparrow}^{\dagger}c_{i\downarrow}\rangle\\
	 *          \langle c_{i\downarrow}^{\dagger}c_{u\uparrow}\rangle       & \langle c_{i\downarrow}^{\dagger}c_{i\downarrow}\rangle
	 *      \end{array}\right].
	 *  \f] */
	virtual Property::Magnetization calculateMagnetization(
		Index pattern,
		Index ranges
	);

	/** Calculate the Magnetization. This function should be overriden by
	 *  those deriving classes that provide support for calculating the
	 *  magnetization. By default the PropertyExtractor prints an error
	 *  message that the given property is not supported.
	 *
	 *  @param patterns A list of patterns that will be matched against the
	 *  @link Index Indices @endlink in the Model to determine which @link
	 *  Index Indices @endlink for which to calculate the Magnetization.
	 *
	 *  @return A Property::Magnetization for the @link Index Indices
	 *  @endlink that match the patterns. */
	virtual Property::Magnetization calculateMagnetization(
		std::vector<Index> patterns
	);

	/** Calculate the local density of states. This function should be
	 *  overriden by those deriving classes that provide support for
	 *  calculating the local density of states. By default the
	 *  PropertyExtractor prints an error message that the given property
	 *  is not supported.
	 *
	 *  @param pattern Specifies the index pattern for which to calculate
	 *  the LDOS. For example, assume that the index scheme is
	 *  {x, y, z, spin}. {ID_X, 5, 10, IDX_SUM_ALL} will calculate the
	 *  LDOS for each x along (y,z)=(5,10) by summing over spin. Similarly
	 *  {ID_X, 5, IDX_Y, IDX_SUM_ALL} will return a two dimensional LDOS
	 *  for all x and z and y = 5. Note that IDX_X, IDX_Y, and IDX_Z refers
	 *  to the first, second, and third index used by the routine to create
	 *  a one-, two-, or three-dimensional output, rather than being tied
	 *  to the x, y, and z used as physical subindices.
	 *
	 *  @param ranges Speifies the number of elements for each subindex. Is
	 *  ignored for indices specified with positive integers in the
	 *  pattern, but is used to loop from 0 to the value in ranges for
	 *  IDX_X, IDX_Y, IDX_Z, and IDX_SUM_ALL. Appropriate ranges
	 *  corresponding to the two pattern examples above are
	 *  {SIZE_X, 1, 1, NUM_SPINS} and {SIZE_X, 1, SIZE_Z, NUM_SPINS},
	 *  respectively.
	 *
	 *  @return A density array with size equal to the number of points
	 *  included by specified patter-range combination. */
	virtual Property::LDOS calculateLDOS(Index pattern, Index ranges);

	/** Calculate the local density of states. This function should be
	 *  overriden by those deriving classes that provide support for
	 *  calculating the local density of states. By default the
	 *  PropertyExtractor prints an error message that the given property
	 *  is not supported.
	 *
	 *  @param patterns A list of patterns that will be matched against the
	 *  @link Index Indices @endlink in the Model to determine which @link
	 *  Index Indices @endlink for which to calculate the local density of
	 *  states.
	 *
	 *  @return A Property::LDOS for the @link Index Indices @endlink that
	 *  match the patterns. */
	virtual Property::LDOS calculateLDOS(
		std::vector<Index> patterns
	);

	/** Calculate the spin-polarized local density of states. This function
	 *  should be overriden by those deriving classes that provide support
	 *  for calculating the spin-polarized local density of states. By
	 *  default the PropertyExtractor prints an error message that the
	 *  given property is not supported.
	 *
	 *  @param pattern Specifies the index pattern for which to calculate
	 *  the spin-polarized LDOS. For example, assume that the index scheme
	 *  is {x, y, z, spin}. {ID_X, 5, 10, IDX_SPIN} will calculate the
	 *  spin-polarized LDOS for each x along (y,z)=(5,10). Similarly
	 *  {ID_X, 5, IDX_Y, IDX_SPIN} will return a two dimensional
	 *  spin-polarized LDOS for all x and z and y = 5. Note that IDX_X,
	 *  IDX_Y, and IDX_Z refers to the first, second, and third index used
	 *  by the routine to create a one-, two-, or three-dimensional output,
	 *  rather than being tied to the x, y, and z used as physical
	 *  subindices.
	 *
	 *  @param ranges Speifies the number of elements for each subindex. Is
	 *  ignored for indices specified with positive integers in the
	 *  pattern, but is used to loop from 0 to the value in ranges for
	 *  IDX_X, IDX_Y, IDX_Z, and IDX_SUM_ALL. Appropriate ranges
	 *  corresponding to the two pattern examples above are
	 *  {SIZE_X, 1, 1, NUM_SPINS} and {SIZE_X, 1, SIZE_Z, NUM_SPINS},
	 *  respectively.
	 *
	 *  @return A spin-polarized LDOS array with size equal to four times
	 *  the number of points included by specified patter-range
	 *  combination.
	 *  The four entries are
	 *  \f[
	 *      \left[\begin{array}{cc}
	 *          0   & 1\\
	 *          2   & 3
	 *      \end{array}\right] =
	 *      \left[\begin{array}{cc}
	 *          \rho_{i\uparrow i\uparrow}(E)       & \rho_{i\uparrow i\downarrow}(E)\\
	 *          \rho_{i\downarrow i\uparrow}(E)     & \rho_{i\downarrow i\downarrow}(E)\\
	 *      \end{array}\right],
	 *  \f]
	 *  where
	 *  \f[
	 *      \rho_{i\sigma i\sigma'}(E) = \sum_{E_n}\langle\Psi_n|c_{i\sigma}^{\dagger}c_{i\sigma'}|\Psi_n\rangle\delta(E - E_n) .
	 *  \f] */
	virtual Property::SpinPolarizedLDOS calculateSpinPolarizedLDOS(
		Index pattern,
		Index ranges
	);

	/** Calculate the spin-polarized local density of states. This function
	 *  should be overriden by those deriving classes that provide support
	 *  for calculating the spin-polarized local density of states. By
	 *  default the PropertyExtractor prints an error message that the
	 *  given property is not supported.
	 *
	 *  @param patterns A list of patterns that will be matched against the
	 *  @link Index Indices @endlink in the Model to determine which @link
	 *  Index Indices @endlink for which to calculate the spin-polarized
	 *  local density of states.
	 *
	 *  @return A Property::SpinPolarizedLDOS for the @link Index Indices
	 *  @endlink that match the patterns. */
	virtual Property::SpinPolarizedLDOS calculateSpinPolarizedLDOS(
		std::vector<Index> patterns
	);

	/** Calculate the expectation value \f$\langle
	 *  c_{to}^{\dagger}c_{from}\f$. This function should be overriden by
	 *  those deriving classes that provide support for calculating the
	 *  expecation value. By default the PropertyExtractor prints an error
	 *  message that the given property is not supported.
	 *
	 *  @param to The Index on the left operator.
	 *  @param from The index on the right operator.
	 *
	 *  @return The expectation value \f$\langle
	 *  c_{to}^{\dagger}c_{from}\f$. */
	virtual std::complex<double> calculateExpectationValue(
		Index to,
		Index from
	);

	/** Calculate the density of states. This function should be overriden
	 *  by those deriving classes that provide support for calculating the
	 *  density of states. By default the PropertyExtractor prints an error
	 *  message that the given property is not supported.
	 *
	 *  @return A Property::DOS containing the density of states. */
	virtual Property::DOS calculateDOS();

	/** Sample the DOS by averaging over the LDOS for multiple random @link
	 *  Index Indices@endlink. The resulting DOS is normalized to integrate
	 *  to the total number of states in the sample space covered by the
	 *  patterns @link Index Indices@endlink.
	 *
	 *  @param numSamples The number of samples to use.
	 *  @param patterns A list of patterns to randomize over. If the list
	 *  is empty (default value), all @link Index Indices@endlink are
	 *  randomized over.
	 *
	 *  @param seed Seed to use for randomization.
	 *
	 *  @return A Property::DOS containing the sampled DOS. */
	virtual Property::DOS sampleDOS(
		unsigned int numSamples,
		const std::vector<Index> &patterns = {},
		unsigned int seed = time(nullptr)
	);

	/** Calculate the entropy. This function should be overriden by those
	 *  deriving classes that provide support for calculating the entropy.
	 *  By default the PropertyExtractor prints an error message that the
	 *  given property is not supported.
	 *
	 *  @return The entropy. */
	virtual double calculateEntropy();
protected:
	/** Energy type. */
	enum class EnergyType{Real, Matsubara};

	/** Get the energy type.
	 *
	 *  @return The EnergyType for the energy window. */
	EnergyType getEnergyType() const;

	/** Get the energy resolution.
	 *
	 *  @return The energy resolution for the energy window. */
	int getEnergyResolution() const;

	/** Get lower bound for the energy window.
	 *
	 *  @return The lower bound for the energy window. */
	double getLowerBound() const;

	/** Get the upper bound for the energy window.
	 *
	 *  @return The upper bound for the energy window. */
	double getUpperBound() const;

	/** Get the lower Fermionic Matsubara energy index.
	 *
	 *  @return The lower Fermionic Matsubara energy index. */
	int getLowerFermionicMatsubaraEnergyIndex() const;

	/** Get the upper Fermionic Matsubara energy index.
	 *
	 *  @return The upper Fermionic Matsubara energy index. */
	int getUpperFermionicMatsubaraEnergyIndex() const;

	/** Get the lower Bosonic Matsubara energy index.
	 *
	 *  @return The lower Bosonic Matsubara energy index. */
	int getLowerBosonicMatsubaraEnergyIndex() const;

	/** Get the upper Bosonic Matsubara energy index.
	 *
	 *  @return The upper Bosonic Matsubara energy index. */
	int getUpperBosonicMatsubaraEnergyIndex() const;

	/** Info class that is able to pass the most common parameters between
	 *  calculate-functions and the corresponding callbacks. Calculations
	 *  requiring more advanced features should override this class and
	 *  perform a type cast in the callbacks. */
	class Information{
	public:
		/** Constructs a
		 *  PropertyExtractor::PropertyExtractor::Information. */
		Information();

		/** Set the subindex for the spin index.
		 *
		 *  @param spinIndex The subindex that is the spin index. */
		void setSpinIndex(int spinIndex);

		/** Get the subindex for the spin index.
		 *
		 *  @return The spin index. Returns -1 if no spin-index has
		 *  been set. */
		int getSpinIndex() const;
	private:
		/** Spin index. */
		int spinIndex;
	};

	/*** Get the nergy infinitesimal.
	 *
	 *  @return The energy infinitesimal. */
	double getEnergyInfinitesimal() const;

	/** Loops over range indices and calls the given callback function to
	 *  calculate the correct quantity. The function recursively calls
	 *  itself replacing any IDX_SUM_ALL, IDX_X, IDX_Y, and IDX_Z
	 *  specifiers by actual subindices in the range [0, ranges[s]), where
	 *  s is the subindex at which the specifier appears. For example, the
	 *  pattern ranges pair {IDX_SUM_ALL, 2, IDX_X} and {2, 1, 3} will
	 *  result in the callback being called for {0, 2, 0}, {0, 2, 1}, {0,
	 *  2, 2}, {1, 2, 0}, {1, 2, 1}, and {1, 2, 2}. The first and fourth,
	 *  second and fifth, and third and sixth Index will further be passed
	 *  to the callback with the same memory offset since their result
	 *  should be summed.
	 *
	 *  The memory offset is further calculated by traversing the
	 *  subindices of the apttern from right to left and multiplying the
	 *  current offset multiplier by the number of indices in the range
	 *  size for the given subindex. This results in an offset that places
	 *  the elements in consequtive order in increasing order of the Index
	 *  order. Where an Index is considered to come before another Index if
	 *  the first subindex to differ between two @link Index Indices
	 *  @endlink from the left is smaller than the other Index.
	 *
	 *  @param callback A callback function that is called to perform the
	 *  actual calculation for a given Index.
	 *
	 *  @param property Reference to Property where the result is to be
	 *  stored.
	 *
	 *  @param pattern An Index specifying the pattern for which to perform
	 *  the calculation.
	 *
	 *  @param ranges The upper limit (exclusive) for which subindices with
	 *  wildcard specifiers will be replaced. The lower limit is 0.
	 *
	 *  @param currentOffset The memory offset calculated for the given
	 *  pattern Index. Should be zero for the initial call to the function.
	 *
	 *  @param offsetMultiplier Number indicating the block size associated
	 *  with increasing the current subindex by one. Should be equal to the
	 *  number of data elements per Index for the initial call to the
	 *  function.
	 *
	 *  @param information Allows for custom information to be passed
	 *  between the calculate-functions and the correpsonding callbacks. */
	template<typename DataType>
	void calculate(
		void (*callback)(
			PropertyExtractor *cb_this,
			Property::Property &property,
			const Index &index,
			int offset,
			Information &information
		),
		Property::AbstractProperty<DataType> &property,
		Index pattern,
		const Index &ranges,
		int currentOffset,
		int offsetMultiplier,
		Information &information
	);

	/** Loops over the indices satisfying the specified patterns and calls
	 *  the appropriate callback function to calculate the correct
	 *  quantity.
	 *
	 *  @param callback A callback function that is called to perform the
	 *  actual calculation for a given Index.
	 *
	 *  @param allIndices An IndexTree containing all the Indices for which
	 *  the callback should be called.
	 *
	 *  @param memoryLayout The memory layout used for the Property.
	 *  @param abstractProperty The Property that is being calculated.
	 *  @param information Allows for custom information to be passed
	 *  between the calculate-functions and the correpsonding callbacks. */
	template<typename DataType>
	void calculate(
		void (*callback)(
			PropertyExtractor *cb_this,
			Property::Property &property,
			const Index &index,
			int offset,
			Information &information
		),
		const IndexTree &allIndices,
		const IndexTree &memoryLayout,
		Property::AbstractProperty<DataType> &abstractProperty,
		Information &information
	);

	/** Ensure that range indices are on compliant format. I.e., sets the
	 *  range to  one for indices with non-negative pattern value.
	 *
	 *  @param pattern The pattern.
	 *  @param ranges The ranges that will have its subindices set to one
	 *  for every pattern subindex that is non negative. */
	void ensureCompliantRanges(const Index &pattern, Index &ranges);

	/** Extract ranges for loop indices. The subindices with IDX_X, IDX_Y
	 *  and IDX_Z are identified and counted and an array of the same size
	 *  as the number of loop indices is created and filled with the ranges
	 *  for the corrsponding loop subindices.
	 *
	 *  @param pattern A pattern.
	 *  @param ranges The ranges for the given pattern.
	 *  @param loopDimensions Pointer to int that will hold the number of
	 *  loop dimensions after the call has completed.
	 *
	 *  @return A vector that contains the ranges for the loop subindices.
	 */
	std::vector<int> getLoopRanges(
		const Index &pattern,
		const Index &ranges
	);

	/** Generate an IndexTree containing all the @link Index Indices
	 *  @endlink in the HoppingAmplitudeSet that matches the given
	 *  patterns. Before being added to the IndexTree, the @link Index
	 *  Indices @endlink may be modified to replace subindices by their
	 *  corresponding pattern value. I.e. A summation or spin subindex may
	 *  still be labeld such in the IndexTree depending on the flags that
	 *  are passed to the function.
	 *
	 *  The pattern can also be a compund Index consisting of two Indices,
	 *  in which case the pattern matching is applied to each component
	 *  Index separately.
	 *
	 *  @param patterns List of patterns to match against.
	 *  @param The HoppingAmplitudeSet cntaining all the @link Index
	 *  Indices @endlink that will be matched against the patterns.
	 *
	 *  @param keepSummationWildcards If true, summation wildcards in the
	 *  pattern will be preserved in the IndexTree.
	 *
	 *  @param keepSpinWildcards If true, spin wildcards in the pattern
	 *  will be preserved in the IndexTree. */
	IndexTree generateIndexTree(
		std::vector<Index> patterns,
		const HoppingAmplitudeSet &hoppingAmplitudeSet,
		bool keepSummationWildcards,
		bool keepSpinWildcards
	);
private:
	/** Energy type used for energy dependent quantities. */
	EnergyType energyType;

	/** Default energy resolution. */
	static constexpr int ENERGY_RESOLUTION = 1000;

	/** Energy resolution used for energy dependent quantities when the
	 *  energy type is EnergyType::Real. */
	int energyResolution;

	/** Default lower bound. */
	static constexpr double LOWER_BOUND = -1.;

	/** Lower bound used for energy dependent quantities when the energy
	 *  type is EnergyType::Real. */
	double lowerBound;

	/** Default upper bound. */
	static constexpr double UPPER_BOUND = 1.;

	/** Upper bound used for energy dependent quantities when the energy
	 *  type is EnergyType::Real. */
	double upperBound;

	/** Default lower Fermionic Matsubara energy index. */
	static constexpr int LOWER_FERMIONIC_MATSUBARA_ENERGY_INDEX = -1;

	/** Lower Fermionic Matsubara energy index used for Fermionic energies
	 *  when the energy type is EnergyType::Matsubara. */
	int lowerFermionicMatsubaraEnergyIndex;

	/** Default upper Fermionic Matsubara energy index. */
	static constexpr int UPPER_FERMIONIC_MATSUBARA_ENERGY_INDEX = 1;

	/** Upper Fermionic Matsubara energy index used for Fermionic energies
	 *  when the energy type is EnergyType::Matsubara. */
	int upperFermionicMatsubaraEnergyIndex;

	/** Default lower Bosonic Matsubara energy index. */
	static constexpr int LOWER_BOSONIC_MATSUBARA_ENERGY_INDEX = 0;

	/** Lower Bosonic Matsubara energy index used for Bosonic energies when
	 *   the energy type is EnergyType::Matsubara. */
	int lowerBosonicMatsubaraEnergyIndex;

	/** Default upper Bosonic Matsubara energy index. */
	static constexpr int UPPER_BOSONIC_MATSUBARA_ENERGY_INDEX = 0;

	/** Upper Bosonic Matsubara energy index used for Bosonic energies
	 *  when the energy type is EnergyType::Matsubara. */
	int upperBosonicMatsubaraEnergyIndex;

	/** Default energy infinitesimal. */
	static constexpr double ENERGY_INFINITESIMAL = 1e-3;

	/** The nergy infinitesimal \f$\delta\f$ that for example can be used
	 *  in the denominator of the Green's function as \f$i\delta\f$. */
	double energyInfinitesimal;

	/** Should return the solver the PropertyExtractor is using. */
	virtual const Solver::Solver& getSolver() const;
};

inline PropertyExtractor::EnergyType PropertyExtractor::getEnergyType() const{
	return energyType;
}

inline int PropertyExtractor::getEnergyResolution() const{
	MyTBTKAssert(
		energyType == EnergyType::Real,
		"PropertyExtractor::PropertyExtractor::getEnergyResolution()",
		"The energy resolution cannot be accessed when the energy type"
		<< " is Matsubara.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with three arguments if real energies are wanted for the"
		<< " PropertyExtractor."
	);

	return energyResolution;
}

inline double PropertyExtractor::getLowerBound() const{
	MyTBTKAssert(
		energyType == EnergyType::Real,
		"PropertyExtractor::PropertyExtractor::getLowerBound()",
		"The lower bound cannot be accessed when the energy type is"
		<< " Matsubara.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with three arguments if real energies are wanted for the"
		<< " PropertyExtractor."
	);

	return lowerBound;
}

inline double PropertyExtractor::getUpperBound() const{
	MyTBTKAssert(
		energyType == EnergyType::Real,
		"PropertyExtractor::PropertyExtractor::getUpperBound()",
		"The upper bound cannot be accessed when the energy type is"
		<< " Matsubara.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with three arguments if real energies are wanted for the"
		<< " PropertyExtractor."
	);

	return upperBound;
}

inline int PropertyExtractor::getLowerFermionicMatsubaraEnergyIndex() const{
	MyTBTKAssert(
		energyType == EnergyType::Matsubara,
		"PropertyExtractor::PropertyExtractor::getLowerFermionicMatsubaraEnergyIndex()",
		"The lower Fermionic Matsubara energy index cannot be accessed"
		<< " when the energy type is real.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with four arguments if Matsubara energies are wanted for"
		<< " the PropertyExtractor."
	);

	return lowerFermionicMatsubaraEnergyIndex;
}

inline int PropertyExtractor::getUpperFermionicMatsubaraEnergyIndex() const{
	MyTBTKAssert(
		energyType == EnergyType::Matsubara,
		"PropertyExtractor::PropertyExtractor::getUpperFermionicMatsubaraEnergyIndex()",
		"The upper Fermionic Matsubara energy index cannot be accessed"
		<< " when the energy type is real.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with four arguments if Matsubara energies are wanted for"
		<< " the PropertyExtractor."
	);

	return upperFermionicMatsubaraEnergyIndex;
}

inline int PropertyExtractor::getLowerBosonicMatsubaraEnergyIndex() const{
	MyTBTKAssert(
		energyType == EnergyType::Matsubara,
		"PropertyExtractor::PropertyExtractor::getLowerBosonicMatsubaraEnergyIndex()",
		"The lower Bosonic Matsubara energy index cannot be accessed"
		<< " when the energy type is real.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with four arguments if Matsubara energies are wanted for"
		<< " the PropertyExtractor."
	);

	return lowerBosonicMatsubaraEnergyIndex;
}

inline int PropertyExtractor::getUpperBosonicMatsubaraEnergyIndex() const{
	MyTBTKAssert(
		energyType == EnergyType::Matsubara,
		"PropertyExtractor::PropertyExtractor::getUpperBosonicMatsubaraEnergyIndex()",
		"The upper Bosonic Matsubara energy index cannot be accessed"
		<< " when the energy type is real.",
		"Use PropertyExtractor::PropertyExtractor::setEnergyWindow()"
		<< " with four arguments if Matsubara energies are wanted for"
		<< " the PropertyExtractor."
	);

	return upperBosonicMatsubaraEnergyIndex;
}

inline double PropertyExtractor::getEnergyInfinitesimal() const{
	return energyInfinitesimal;
}

template<typename DataType>
void PropertyExtractor::calculate(
	void (*callback)(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	),
	Property::AbstractProperty<DataType> &property,
	Index pattern,
	const Index &ranges,
	int currentOffset,
	int offsetMultiplier,
	Information &information
){
	//Find the next specifier index.
	int currentSubindex = pattern.getSize()-1;
	for(; currentSubindex >= 0; currentSubindex--){
		if(pattern.at(currentSubindex) < 0)
			break;
	}

	if(currentSubindex == -1){
		//No further specifier index found. Call the callback.
		callback(this, property, pattern, currentOffset, information);
	}
	else{
		//Ensure that the specifier is valid for the Ranges format.
		MyTBTKAssert(
			pattern.at(currentSubindex).isSummationIndex()
			|| pattern.at(currentSubindex).isRangeIndex(),
			"PropertyExtractor::calculate()",
			"Invalid specifier found at subindex "
			<< currentSubindex << ".",
			"Did you mean IDX_SUM_ALL, IDX_X, IDX_Y, or IDX_Z?"
		);

		//Multiply the memory offset for non summation indices.
		int nextOffsetMultiplier = offsetMultiplier;
		if(!pattern.at(currentSubindex).isSummationIndex())
			nextOffsetMultiplier *= ranges.at(currentSubindex);

		//Set flag indicating whether the current subindex is a
		//summation index.
		bool isSumIndex = false;
		if(pattern.at(currentSubindex).isSummationIndex())
			isSumIndex = true;

		//Recurively call the calculate function with the specifier at
		//the current subindex replaced by each subindex value in the
		//corresponding range.
		for(int n = 0; n < ranges.at(currentSubindex); n++){
			pattern.at(currentSubindex) = n;
			calculate(
				callback,
				property,
				pattern,
				ranges,
				currentOffset,
				nextOffsetMultiplier,
				information
			);
			if(!isSumIndex)
				currentOffset += offsetMultiplier;
		}
	}
}

template<typename DataType>
void PropertyExtractor::calculate(
	void (*callback)(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	),
	const IndexTree &allIndices,
	const IndexTree &memoryLayout,
	Property::AbstractProperty<DataType> &abstractProperty,
	Information &information
){
	for(
		IndexTree::ConstIterator iterator = allIndices.cbegin();
		iterator != allIndices.end();
		++iterator
	){
		const Index &index = *iterator;
		std::vector<unsigned int> spinIndices
			= memoryLayout.getSubindicesMatching(
				IDX_SPIN,
				index,
				IndexTree::SearchMode::MatchWildcards
			);
		if(spinIndices.size() != 0){
			MyTBTKAssert(
				spinIndices.size() == 1,
				"PropertyExtractor::calculate()",
				"Several spin indeces found.",
				"Use IDX_SPIN at most once per pattern to"
				<< " indicate spin index."
			);
			information.setSpinIndex(spinIndices[0]);
		}

		callback(
			this,
			abstractProperty,
			index,
			abstractProperty.getOffset(index),
			information
		);
	}
}

inline void PropertyExtractor::Information::setSpinIndex(int spinIndex){
	this->spinIndex = spinIndex;
}

inline int PropertyExtractor::Information::getSpinIndex() const{
	return spinIndex;
}

inline const Solver::Solver& PropertyExtractor::getSolver() const{
	MyTBTKExit(
		"PropertyExtractor::getSolver()",
		"Missing implementation of PropertyExtractor::getSolver() in"
		" deriving class.",
		""
	);
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
