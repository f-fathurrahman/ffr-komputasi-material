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

/// @cond MyTBTK_FULL_DOCUMENTATION
/** @package MyTBTKcalc
 *  @file LadderOperator.h
 *  @brief Ladder operator.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_LADDER_OPERATOR
#define COM_DAFER45_MyTBTK_LADDER_OPERATOR

#include "MyTBTK/FockState.h"
#include "MyTBTK/HoppingAmplitudeSet.h"
#include "MyTBTK/Statistics.h"

namespace MyTBTK{

template<typename BIT_REGISTER>
class LadderOperator{
public:
	/** Operator type. */
	enum class Type {Creation, Annihilation};

	/** Constructor. */
	LadderOperator();

	/** Constructor. */
	LadderOperator(
		Type type,
		Statistics statistics,
		const HoppingAmplitudeSet *hoppingAmplitudeSet,
		unsigned int state,
		unsigned int numBitsPerState,
		unsigned int maxOccupation,
		const FockState<BIT_REGISTER> &templateState,
		const BIT_REGISTER &fermionMask
	);

	/** Destructor. */
	~LadderOperator();

	/** Get type. */
	Type getType() const;

	/** Get physical Index of associated single-particle state. */
	const Index getPhysicalIndex() const;

	/** Get associated single-particle state. */
	unsigned int getState() const;

	/** Get number of particles in the associated state. */
	unsigned int getNumParticles(
		const FockState<BIT_REGISTER> &fockState
	) const;

	/** Multiplication operator. */
	FockState<BIT_REGISTER>& operator*(FockState<BIT_REGISTER> &rhs) const;
private:
	/** Operator type. */
	Type type;

	/** Operator statistics. */
	Statistics statistics;

	/** Pointer to the corresponding HoppingAmplitudeSet. */
	const HoppingAmplitudeSet *hoppingAmplitudeSet;

	/** Single-particle state index. */
	unsigned int state;

	/** State mask. */
	BIT_REGISTER stateMask;

	/** Least significant bit. */
	BIT_REGISTER leastSignificantBit;

	/** Index of least significant bit. */
	unsigned int leastSignificantBitIndex;

	/** State corresponding to maximum number of occupied particles. */
	BIT_REGISTER maxOccupation;

	/** Mask for singeling out those fermions that have a higher bit index
	 *  than the state corresponding to this opperator. */
	BIT_REGISTER moreSignificantFermionMask;
};

template<typename BIT_REGISTER>
LadderOperator<BIT_REGISTER>::LadderOperator(){
}

template<typename BIT_REGISTER>
LadderOperator<BIT_REGISTER>::LadderOperator(
	Type type,
	Statistics statistics,
	const HoppingAmplitudeSet *hoppingAmplitudeSet,
	unsigned int state,
	unsigned int numBitsPerState,
	unsigned int maxOccupation,
	const FockState<BIT_REGISTER> &templateState,
	const BIT_REGISTER &fermionMask
) :
	stateMask(templateState.bitRegister),
	leastSignificantBit(templateState.bitRegister),
	maxOccupation(templateState.bitRegister),
	moreSignificantFermionMask(fermionMask)
{
	this->type = type;
	this->statistics = statistics;
	this->hoppingAmplitudeSet = hoppingAmplitudeSet;
	this->state = state;

	leastSignificantBitIndex = state*numBitsPerState;

	for(unsigned int n = 0; n < stateMask.getNumBits(); n++){
		if(n >= leastSignificantBitIndex && n < leastSignificantBitIndex + numBitsPerState)
			stateMask.setBit(n, 1);
		else
			stateMask.setBit(n, 0);

		if(n == leastSignificantBitIndex)
			leastSignificantBit.setBit(n, 1);
		else
			leastSignificantBit.setBit(n, 0);
	}

	this->maxOccupation = maxOccupation;
	this->maxOccupation = (this->maxOccupation << leastSignificantBitIndex);

	for(unsigned int n = 0; n < moreSignificantFermionMask.getNumBits(); n++){
		this->moreSignificantFermionMask.setBit(n, false);
		if(leastSignificantBit.getBit(n))
			break;
	}
}

template<typename BIT_REGISTER>
LadderOperator<BIT_REGISTER>::~LadderOperator(){
}

template<typename BIT_REGISTER>
typename LadderOperator<BIT_REGISTER>::Type LadderOperator<BIT_REGISTER>::getType() const{
	return type;
}

template<typename BIT_REGISTER>
const Index LadderOperator<BIT_REGISTER>::getPhysicalIndex() const{
	return hoppingAmplitudeSet->getPhysicalIndex(state);
}

template<typename BIT_REGISTER>
unsigned int LadderOperator<BIT_REGISTER>::getState() const{
	return state;
}

template<typename BIT_REGISTER>
unsigned int LadderOperator<BIT_REGISTER>::getNumParticles(
	const FockState<BIT_REGISTER> &fockState
) const{
	return ((fockState.getBitRegister() & stateMask) >> leastSignificantBitIndex).toUnsignedInt();
}

template<typename BIT_REGISTER>
FockState<BIT_REGISTER>& LadderOperator<BIT_REGISTER>::operator*(
	FockState<BIT_REGISTER> &rhs
) const{
	switch(type){
	case Type::Creation:
		if((rhs.bitRegister & stateMask) == maxOccupation){
			rhs.bitRegister.setMostSignificantBit();
			break;
		}
		rhs.bitRegister += leastSignificantBit;
		break;
	case Type::Annihilation:
		if(!(rhs.bitRegister & stateMask).toBool()){
			rhs.bitRegister.setMostSignificantBit();
			break;
		}
		rhs.bitRegister -= leastSignificantBit;
		break;
	default:
		MyTBTKExit(
			"LadderOperator<BIT_REGISTER>::operator*()",
			"This should never happen.",
			"Contact the developer."
		);
	}

	switch(statistics){
	case Statistics::FermiDirac:
		rhs.prefactor *= pow(-1, (rhs.bitRegister & moreSignificantFermionMask).getNumOneBits());
		break;
	case Statistics::BoseEinstein:
		break;
	default:
		MyTBTKExit(
			"LadderOperator<BIT_REGISTER>::operator*()",
			"This should never happen.",
			"Contact the developer."
		);
	}

	return rhs;
}

};	//End of namespace MyTBTK

#endif
/// @endcond
