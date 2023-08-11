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
 *  @file SelfEnergy.h
 *  @brief Extracts physical properties from the
 *  Solver::SelfEnergy.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_SELF_ENERGY
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_SELF_ENERGY

#include "MyTBTK/Solver/SelfEnergy.h"
#include "MyTBTK/Property/SelfEnergy.h"
#include "MyTBTK/PropertyExtractor/PropertyExtractor.h"

#include <complex>
//#include <initializer_list>

namespace MyTBTK{
namespace PropertyExtractor{

/** The PropertyExtractor::SelfEnergy extracts the SelfEnergy from
 *  Solver::SelfEnergy. */
class SelfEnergy : public PropertyExtractor{
public:
	/** Constructs a PropertyExtractor::SelfEnergy.
	 *
	 *  @param solver The Solver to use. */
	SelfEnergy(Solver::SelfEnergy &solver);

	/** Calculates the Susceptibility. */
	virtual Property::SelfEnergy calculateSelfEnergy(
		std::vector<Index> patterns
	);
private:
	/** Calback for callculating the self-energy. */
	static void calculateSelfEnergyCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Solver::SelfEnergy to work on. */
	Solver::SelfEnergy *solver;

	/** Energies. */
	std::vector<std::complex<double>> energies;

	/** Implements PropertyExtractor::getSolver(). */
	virtual const Solver::Solver& getSolver() const;
};

inline const Solver::Solver& SelfEnergy::getSolver() const{
	return *solver;
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
/// @endcond
