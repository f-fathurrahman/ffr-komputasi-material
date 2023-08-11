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
 *  @file LookupTableMap.h
 *  @brief LookupTableMap.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_LOOKUP_TABLE_MAP
#define COM_DAFER45_MyTBTK_LOOKUP_TABLE_MAP

#include "MyTBTK/FockStateMap/FockStateMap.h"
#include "MyTBTK/BitRegister.h"
#include "MyTBTK/ExtensiveBitRegister.h"

namespace MyTBTK{
namespace FockStateMap{

template<typename BIT_REGISTER>
class LookupTableMap : public FockStateMap<BIT_REGISTER>{
public:
	/** Constructor. */
	LookupTableMap(unsigned int exponentialDimension);

	/** Destructor. */
	virtual ~LookupTableMap();

	/** Get many-body Hilbert space size. */
	virtual unsigned int getBasisSize() const;

	/** Get many-body Hilbert space index for corresponding FockState. */
	virtual unsigned int getBasisIndex(
		const FockState<BIT_REGISTER> &fockState
	) const;

	/** Get FockState for corresponding many-body Hilbert space index. */
	virtual FockState<BIT_REGISTER> getFockState(unsigned int index) const;

	/** Add state. */
	void addState(const FockState<BIT_REGISTER> &fockState);
private:
	/** List of FockStates. */
	std::vector<FockState<BIT_REGISTER>> states;
};

template<typename BIT_REGISTER>
LookupTableMap<BIT_REGISTER>::LookupTableMap(
	unsigned int exponentialDimension
) :
	FockStateMap<BIT_REGISTER>(exponentialDimension)
{
}

template<typename BIT_REGISTER>
LookupTableMap<BIT_REGISTER>::~LookupTableMap(){
}

template<typename BIT_REGISTER>
unsigned int LookupTableMap<BIT_REGISTER>::getBasisSize() const{
	return states.size();
}

template<typename BIT_REGISTER>
unsigned int LookupTableMap<BIT_REGISTER>::getBasisIndex(
	const FockState<BIT_REGISTER> &fockState
) const{
	unsigned int min = 0;
	unsigned int max = states.size()-1;
	while(min <= max){
		unsigned int currentState = (min+max)/2;
		if(fockState.getBitRegister() > states.at(currentState).getBitRegister())
			min = currentState + 1;
		else if(fockState.getBitRegister() < states.at(currentState).getBitRegister())
			max = currentState - 1;
		else if(fockState.getBitRegister() == states.at(currentState).getBitRegister())
			return currentState;
	}
	MyTBTKExit(
		"LookupTableFockStateMap<BIT_REGISTER>::getBasisIndex()",
		"FockState not found.",
		""
	);
}

template<typename BIT_REGISTER>
FockState<BIT_REGISTER> LookupTableMap<BIT_REGISTER>::getFockState(
	unsigned int index
) const{
	return states.at(index);
}

template<typename BIT_REGISTER>
void LookupTableMap<BIT_REGISTER>::addState(
	const FockState<BIT_REGISTER> &fockState
){
	if(
		states.size() == 0
		|| states.back().getBitRegister() < fockState.getBitRegister()
	){
		states.push_back(fockState);
	}
	else{
		unsigned int min = 0;
		unsigned int max = states.size()-1;
		while(min <= max){
			unsigned int currentState = (min+max)/2;
			if(
				fockState.getBitRegister()
				> states.at(currentState).getBitRegister()
			){
				min = currentState + 1;
			}
			else if(
				fockState.getBitRegister()
				< states.at(currentState).getBitRegister()
			){
				max = currentState - 1;
			}

			if(min >= max){
				if(
					fockState.getBitRegister()
						< states.at(
							currentState
						).getBitRegister()
				){
					states.insert(
						states.begin() + currentState,
						fockState
					);
				}
				else{
					states.insert(
						states.begin() + currentState
							+ 1,
						fockState
					);
				}
				break;
			}
		}
	}
}

};	//End of namespace FockStateMap
};	//End of namespace MyTBTK

#endif
/// @endcond
