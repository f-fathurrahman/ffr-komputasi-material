/* Copyright 2017 Kristofer Björnson
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
 *  @file BlockDiagonalizer.h
 *  @brief Extracts physical properties from the Solver::BlockDiagonalizer.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_BLOCK_DIAGONALIZER
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_BLOCK_DIAGONALIZER

#include "MyTBTK/Solver/BlockDiagonalizer.h"
#include "MyTBTK/Property/DOS.h"
#include "MyTBTK/Property/Density.h"
#include "MyTBTK/Property/EigenValues.h"
#include "MyTBTK/Property/GreensFunction.h"
#include "MyTBTK/Property/LDOS.h"
#include "MyTBTK/Property/Magnetization.h"
#include "MyTBTK/Property/SpinPolarizedLDOS.h"
#include "MyTBTK/Property/WaveFunctions.h"
#include "MyTBTK/PropertyExtractor/PropertyExtractor.h"

#include <complex>
//#include <initializer_list>

namespace MyTBTK{
namespace PropertyExtractor{

/** @brief Extracts physical properties from the Solver::BlockDiagonalizer.
 *
 *  The PropertyExtractor::Diagonalizer extracts @link
 *  Property::AbstractProperty@endlink from the Solver::BlockDiagonalizer.
 *
 *  # Example
 *  \snippet PropertyExtractor/BlockDiagonalizer.cpp BlockDiagonalizer
 *  ## Output
 *  \snippet output/PropertyExtractor/BlockDiagonalizer.txt BlockDiagonalizer
 *  \image html output/PropertyExtractor/BlockDiagonalizer/figures/PropertyExtractorBlockDiagonalizerDOS.png
 *  \image html output/PropertyExtractor/BlockDiagonalizer/figures/PropertyExtractorBlockDiagonalizerDensity.png */
class BlockDiagonalizer : public PropertyExtractor{
public:
	/** Constructs a PropertyExtractor::BlockDiagonalizer.
	 *
	 *  @param solver The Solver to use. */
	BlockDiagonalizer(Solver::BlockDiagonalizer &solver);

	/** Get eigenvalues. The eigenvalues are ordered first by block, and
	 *  then in accending order. This means that eigenvalues for blocks
	 *  with smaller @link Index Indices @endlink comes before eigenvalues
	 *  for blocks with larger @link Index Indices @endlink, while inside
	 *  each block the eigenvalues are in accending order.
	 *
	 *  @return A Property::EigenValues containing all the eigenvalues. */
	Property::EigenValues getEigenValues();

	/** Get eigenvalue. The eigenvalues are ordered first by block, and
	 *  then in accending order. This means that eigenvalues for blocks
	 *  with smaller @link Index Indices @endlink comes before eigenvalues
	 *  for blocks with larger @link Index Indices @endlink, while inside
	 *  each block the eigenvalues are in accending order.
	 *
	 *  @param state The state number.
	 *
	 *  @return The eigenvalue for the given state. */
	double getEigenValue(int state) const;

	/** Get eigenvalue. The eigenvalues are ordered in accending order.
	 *
	 *  @param blockIndex The block Index to get the eigenvalue for.
	 *  @param state The intra block state index starting from zero for the
	 *  lowest eigenvalue in the block.
	 *
	 *  @return The eigenvalue for the given block and state. */
	double getEigenValue(const Index &blockIndex, int state) const;

	/** Get amplitude for given eigenvector \f$n\f$ and physical index
	 *  \f$x\f$: \f$\Psi_{n}(x)\f$.
	 *
	 *  @param state Eigenstate number \f$n\f$
	 *  @param index Physical index \f$x\f$.
	 *
	 *  @return The amplitude \f$\Psi_{n}(\mathbf{x})\f$. */
	const std::complex<double> getAmplitude(int state, const Index &index);

	/** Get amplitude for given eigenvector \f$n\f$, block index \f$b\f$,
	 *  and physical Index \f$x\f$: \f$\Psi_{nb}(\mathbf{x})\f$.
	 *
	 *  @param blockIndex Block Index \f$b\f$.
	 *  @param state Eigenstate number \f$n\f$ relative to the given block.
	 *  @param index Physical Index \f$x\f$.
	 *
	 *  @return The amplitude \f$\Psi_{nb}(\mathbf{x})\f$.  */
	const std::complex<double> getAmplitude(
		const Index &blockIndex,
		int state,
		const Index &intraBlockIndex
	) const;

	/** Calculate the wave function on the Custom format. [See
	 *  AbstractProperty for detailed information about the Custom format.
	 *  See PropertyExtractor for detailed information about the patterns
	 *  argument.]
	 *
	 *  @param patterns The pattern to use
	 *  @param states The states to extract the wave function for. Can be
	 *  set to {IDX_ALL} to get all states.
	 *
	 *  @return A WaveFunctions object containing the wave functions values
	 *  for the Indices that satisfies the given patterns and state
	 *  numbers. */
	Property::WaveFunctions calculateWaveFunctions(
		std::vector<Index> patterns,
		std::vector<Subindex> states
	);

	/** Calculate the Green's function on the Custom format. [See
	 *  AbstractProperty for detailed information about the Custom format.
	 *  See PropertyExtractor for detailed information about the patterns
	 *  argument.]
	 *
	 *  @param patterns The pattern to use.
	 *  @param type The Green's function type.
	 *
	 *  @return A GreensFunction for the given patterns. */
	Property::GreensFunction calculateGreensFunction(
		std::vector<Index> patterns,
		Property::GreensFunction::Type type
			= Property::GreensFunction::Type::Retarded
	);

	/** Overrides PropertyExtractor::calculateDOS(). */
	virtual Property::DOS calculateDOS();

	/** Calculate expectation value. */
	virtual std::complex<double> calculateExpectationValue(
		Index to,
		Index from
	);

	/** Overrides PropertyExtractor::calculateDensity(). */
	virtual Property::Density calculateDensity(
		std::vector<Index> patterns
	);

	/** Overrides PropertyExtractor::calculateMagnetization(). */
	virtual Property::Magnetization calculateMagnetization(
		std::vector<Index> patterns
	);

	/** Overrides PropertyExtractor::calculateLDOS(). */
	virtual Property::LDOS calculateLDOS(
		std::vector<Index> patterns
	);

	/** Overrides PropertyExtractor::calculateSpinPolarizedLDOS(). */
	virtual Property::SpinPolarizedLDOS calculateSpinPolarizedLDOS(
		std::vector<Index> patterns
	);

	/** Overrider PropertyExtractor::calculateEntropy(). */
	virtual double calculateEntropy();
private:
	/** Callback for calculating the wave function. Used by
	 *  calculateWaveFunctions. */
	static void calculateWaveFunctionsCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Callback for calculating the Green's function. Used by
	 *  calculateGreensFunction. */
	static void calculateGreensFunctionCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Callback for calculating density. Used by calculateDensity. */
	static void calculateDensityCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Callback for calculating magnetization. Used by calculateMAG. */
	static void calculateMagnetizationCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Calback for callculating local density of states. Used by
	 *  calculateLDOS. */
	static void calculateLDOSCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Callback for calculating spin-polarized local density of states.
	 *  Used by calculateSP_LDOS. */
	static void calculateSP_LDOSCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Solver::Diagonalizer to work on. */
	Solver::BlockDiagonalizer &solver;

	/** Implements PropertyExtractor::getSolver(). */
	virtual const Solver::Solver& getSolver() const;
};

inline double BlockDiagonalizer::getEigenValue(int state) const{
	return solver.getEigenValue(state);
}

inline double BlockDiagonalizer::getEigenValue(
	const Index &blockIndex,
	int state
) const{
	return solver.getEigenValue(blockIndex, state);
}

inline const std::complex<double> BlockDiagonalizer::getAmplitude(
	int state,
	const Index &index
){
	return solver.getAmplitude(state, index);
}

inline const std::complex<double> BlockDiagonalizer::getAmplitude(
	const Index &blockIndex,
	int state,
	const Index &intraBlockIndex
) const{
	return solver.getAmplitude(blockIndex, state, intraBlockIndex);
}

inline const Solver::Solver& BlockDiagonalizer::getSolver() const{
	return solver;
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
