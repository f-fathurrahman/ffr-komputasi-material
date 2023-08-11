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
 *  @file Greens.h
 *  @brief Extracts physical properties from the Solver::Greens.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_GREENS
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_GREENS

#include "MyTBTK/Solver/Greens.h"
#include "MyTBTK/Property/Density.h"
#include "MyTBTK/Property/GreensFunction.h"
#include "MyTBTK/Property/LDOS.h"
#include "MyTBTK/Property/Magnetization.h"
#include "MyTBTK/Property/SpinPolarizedLDOS.h"
#include "MyTBTK/PropertyExtractor/PropertyExtractor.h"

//#include <initializer_list>
#include <iostream>

namespace MyTBTK{
namespace PropertyExtractor{

/** Experimental class for extracting properties from a Solver::Greens. */
class Greens : public PropertyExtractor{
public:
	/** Constructor. */
	Greens(Solver::Greens &cSolver);

	/** Destructor. */
	virtual ~Greens();

	/** Overrides PropertyExtractor::setEnergyWindow(). */
	virtual void setEnergyWindow(
		double lowerBound,
		double upperBound,
		int energyResolution
	);

	/** Overrides PropertyExtractor::calculateDensity(). */
	virtual Property::Density calculateDensity(
		Index pattern,
		Index ranges
	);

	/** Overrides PropertyExtractor::calculateDensity(). */
	virtual Property::Density calculateDensity(
		std::vector<Index> patterns
	);

	/** Overrides PropertyExtractor::calculateLDOS(). */
	virtual Property::LDOS calculateLDOS(
		std::vector<Index> patterns
	);
private:
	/** ChebyshevExpander to work on. */
	Solver::Greens *solver;

	/** Implements PropertyExtractor::getSolver(). */
	virtual const Solver::Solver& getSolver() const;

	/** Callback for calculating the density. Used by calculateDensity. */
	static void calculateDensityCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Callback for calculating the local density of states. Used by
	 *  calculateLDOS. */
	static void calculateLDOSCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);
};

inline const Solver::Solver& Greens::getSolver() const{
	return *solver;
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
/// @endcond
