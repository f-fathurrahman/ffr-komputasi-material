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
 *  @file ElectronFluctuationVertex.h
 *  @brief Extracts physical properties from the
 *  Solver::ElectronFluctuationVertex.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_ELECTRON_FLUCTUATION_VERTEX
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_ELECTRON_FLUCTUATION_VERTEX

#include "MyTBTK/Solver/ElectronFluctuationVertex.h"
#include "MyTBTK/Property/InteractionVertex.h"
#include "MyTBTK/PropertyExtractor/PropertyExtractor.h"

#include <complex>
//#include <initializer_list>

namespace MyTBTK{
namespace PropertyExtractor{

/** The PropertyExtractor::ElectronFluctuationVertex extracts the
 *  InteractionVertex from Solver::ElectronFluctuationVertex. */
class ElectronFluctuationVertex : public PropertyExtractor{
public:
	/** Constructs a PropertyExtractor::ELectronFluctuationVertex.
	 *
	 *  @param solver The Solver to use. */
	ElectronFluctuationVertex(Solver::ElectronFluctuationVertex &solver);

	/** Overrides PropertyExtractor::setEnergyWindow(). */
	virtual void setEnergyWindow(
		double lowerBound,
		double upperBound,
		int resolution
	);

	/** Overrides PropertyExtractor::setEnergyWindow(). */
	virtual void setEnergyWindow(
		int lowerFermionicMatsubaraEnergyIndex,
		int upperFermionicMatsubaraEnergyIndex,
		int lowerBosonicMatsubaraEnergyIndex,
		int upperBosonicEnergyIndex
	);

	/** Calculates the Susceptibility. */
	virtual Property::InteractionVertex calculateInteractionVertex(
		std::vector<Index> patterns
	);
private:
	/** Calback for callculating susceptibility. */
	static void calculateInteractionVertexCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Solver::ElectronFluctautionVertex to work on. */
	Solver::ElectronFluctuationVertex *solver;

	/** Implements PropertyExtractor::getSolver(). */
	virtual const Solver::Solver& getSolver() const;
};

inline const Solver::Solver& ElectronFluctuationVertex::getSolver() const{
	return *solver;
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
/// @endcond
