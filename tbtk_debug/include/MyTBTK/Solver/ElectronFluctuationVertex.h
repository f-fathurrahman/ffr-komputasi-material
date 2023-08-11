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

/// @cond MyTBTK_FULL_DOCUMENTATION
/** @package MyTBTKcalc
 *  @file SelfEnergyCalculator.h
 *  @brief Calculates the self-energy using the RPA approximation.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_SOLVER_ELECTRON_FLUCTUATION_VERTEX
#define COM_DAFER45_MyTBTK_SOLVER_ELECTRON_FLUCTUATION_VERTEX

#include "MyTBTK/BrillouinZone.h"
#include "MyTBTK/IndexedDataTree.h"
#include "MyTBTK/Property/Susceptibility.h"
#include "MyTBTK/MomentumSpaceContext.h"
#include "MyTBTK/Solver/Solver.h"

namespace MyTBTK{
namespace Solver{

class ElectronFluctuationVertex : public Solver, public Communicator{
public:
	/** Constructor. */
	ElectronFluctuationVertex(
		const MomentumSpaceContext &momentumSpaceContext,
		const Property::Susceptibility &susceptibility
	);

	/** Get momentum space context. */
	const MomentumSpaceContext& getMomentumSpaceContext() const;

	/** Get the susceptibility. */
	const Property::Susceptibility& getSusceptibility() const;

	/** Calculate self-energy vertex. */
	std::vector<std::complex<double>> calculateSelfEnergyVertex(
		const Index &index
	);

	/** Set the left interaction. */
	void setLeftInteraction(
		const std::vector<InteractionAmplitude> &leftInteraction
	);

	/** Set the right interaction. */
	void setRightInteraction(
		const std::vector<InteractionAmplitude> &rightInteraction
	);
private:
	/** Momentum space context. */
	const MomentumSpaceContext &momentumSpaceContext;

	/** Susceptibility. */
	const Property::Susceptibility &susceptibility;

	/** Left interaciton. */
	std::vector<InteractionAmplitude> leftInteraction;

	/** Left interaciton. */
	std::vector<InteractionAmplitude> rightInteraction;

	/** Main algorithm for calculating the self-energy vertex.*/
	void calculateSelfEnergyVertexMainAlgorithm(
		std::vector<std::complex<double>> &selfEnergyVertex,
		const Index &index,
		const Property::Susceptibility &susceptibility,
		const std::vector<InteractionAmplitude> &uLeft,
		const std::vector<InteractionAmplitude> &uRight
	);
};

inline const MomentumSpaceContext&
ElectronFluctuationVertex::getMomentumSpaceContext() const{
	return momentumSpaceContext;
}

inline const Property::Susceptibility&
ElectronFluctuationVertex::getSusceptibility() const{
	return susceptibility;
}

inline void ElectronFluctuationVertex::setLeftInteraction(
	const std::vector<InteractionAmplitude> &leftInteraction
){
	this->leftInteraction = leftInteraction;
}

inline void ElectronFluctuationVertex::setRightInteraction(
	const std::vector<InteractionAmplitude> &rightInteraction
){
	this->rightInteraction = rightInteraction;
}

};	//End of namespace Solver
};	//End of namespace MyTBTK

#endif
/// @endcond
