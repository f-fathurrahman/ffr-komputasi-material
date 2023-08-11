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
 *  @file SelfEnergy2.h
 *  @brief Extracts physical properties from the
 *  Solver::SelfEnergy2.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_SELF_ENERGY2
#define COM_DAFER45_MyTBTK_PROPERTY_EXTRACTOR_SELF_ENERGY2

#include "MyTBTK/Solver/SelfEnergy2.h"
#include "MyTBTK/Property/SelfEnergy.h"
#include "MyTBTK/PropertyExtractor/PropertyExtractor.h"

#include <complex>

namespace MyTBTK{
namespace PropertyExtractor{

/** The PropertyExtractor::SelfEnergy extracts the SelfEnergy from
 *  Solver::SelfEnergy2. */
class SelfEnergy2 : public PropertyExtractor{
public:
	/** Constructs a PropertyExtractor::SelfEnergy2.
	 *
	 *  @param solver The Solver to use. */
	SelfEnergy2(Solver::SelfEnergy2 &solver);

	/** Calculates the Susceptibility. */
	virtual Property::SelfEnergy calculateSelfEnergy(
		std::vector<Index> patterns
	);
private:
	/***/
	class SelfEnergyBlockInformation : public Information{
	public:
		/** Constructs a
		 *  PropertyExtractor::SelfEnergy2::SelfEnergyBlockInformation.
		 */
		SelfEnergyBlockInformation();

		/** Set whether the self-energy should be calculated for all
		 *  block indices. */
		void setCalculateSelfEnergyForAllBlocks(
			bool calculateSelfEnergyForAllBlocks
		);

		/** Get whether the self-energy should be calculated for all
		 *  block indices. */
		bool getCalculateSelfEnergyForAllBlocks() const;
	private:
		/** Flag indicating whether the self-energy should be
		 *  calculated for all block indices. */
		bool calculateSelfEnergyForAllBlocks;
	};

	/** Calback for callculating the self-energy. */
	static void calculateSelfEnergyCallback(
		PropertyExtractor *cb_this,
		Property::Property &property,
		const Index &index,
		int offset,
		Information &information
	);

	/** Solver::SelfEnergy2 to work on. */
	Solver::SelfEnergy2 *solver;

	/** Implements PropertyExtractor::getSolver(). */
	virtual const Solver::Solver& getSolver() const;
};

inline void SelfEnergy2::SelfEnergyBlockInformation::setCalculateSelfEnergyForAllBlocks(
	bool calculateSelfEnergyForAllBlocks
){
	this->calculateSelfEnergyForAllBlocks
		= calculateSelfEnergyForAllBlocks;
}

inline bool SelfEnergy2::SelfEnergyBlockInformation::getCalculateSelfEnergyForAllBlocks(
) const{
	return calculateSelfEnergyForAllBlocks;
}

inline const Solver::Solver& SelfEnergy2::getSolver() const{
	return *solver;
}

};	//End of namespace PropertyExtractor
};	//End of namespace MyTBTK

#endif
/// @endcond
