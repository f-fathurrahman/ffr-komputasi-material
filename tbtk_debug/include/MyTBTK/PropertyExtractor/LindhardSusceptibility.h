/* Copyright 2018 Kristofer Björnson
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

/// @cond MyTBTK_FULL_DOCUMENTATION
/** @package MyTBTKcalc
 *  @file DPropertyExtractor.h
 *  @brief Extracts physical properties from the
 *  Solver::LindhardSusceptibility.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_LINDHARD_SUSCEPTIBILITY
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_LINDHARD_SUSCEPTIBILITY

#include "MyTBTK/Solver/LindhardSusceptibility.h"
#include "MyTBTK/Property/Susceptibility.h"
#include "MyTBTK/PropertyExtractor/PropertyExtractor.h"

#include <complex>

namespace MyTBTK{
namespace PropertyExtractor{

/** The PropertyExtractor::LindhardSusceptibility extracts the Susceptibility
 *  from Solver::LindhardSusceptibility. */
class LindhardSusceptibility : public PropertyExtractor{
public:
	/** Constructs a PropertyExtractor::Diagonalizer.
	 *
	 *  @param solver The Solver to use. */
	LindhardSusceptibility(Solver::LindhardSusceptibility &solver);

	/** Calculates the Susceptibility. */
	virtual Property::Susceptibility calculateSusceptibility(
		std::vector<Index> patterns
	);
private:
	/** Calback for callculating susceptibility. */
	static void calculateSusceptibilityCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Solver::Diagonalizer to work on. */
	Solver::LindhardSusceptibility *solver;

	/** Implements PropertyExtractor::getSolver(). */
	virtual const Solver::Solver& getSolver() const;

	/** Energies. */
	std::vector<std::complex<double>> energies;
};

inline const Solver::Solver& LindhardSusceptibility::getSolver() const{
	return *solver;
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
/// @endcond
