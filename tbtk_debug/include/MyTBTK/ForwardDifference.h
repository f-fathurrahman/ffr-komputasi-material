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

/** @package MyTBTKcalc
 *  @file ForwardDifference.h
 *  @brief HoppingAmplitudeList corresponding to a forward difference.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_FORWARD_DIFFERENCE
#define COM_DAFER45_MyTBTK_FORWARD_DIFFERENCE

#include "MyTBTK/HoppingAmplitudeList.h"

namespace MyTBTK{

/** @brief HoppingAmplitudeList corresponding to a forward difference. */
class ForwardDifference : public HoppingAmplitudeList{
public:
	/** Constructs a ForwardDifference.
	 *
	 *  @param subindex Subindex corresponding to the direction of
	 *  differentiation.
	 *
	 *  @param index Index at which the difference is taken.
	 *  @param dx Step length. */
	ForwardDifference(
		unsigned int subindex,
		const Index &index,
		double dx = 1
	);
};

inline ForwardDifference::ForwardDifference(
	unsigned int subindex,
	const Index &index,
	double dx
){
	MyTBTKAssert(
		subindex < index.getSize(),
		"ForwardDifference::ForwardDifference()",
		"Invalid subindex. The subindex '" << subindex << "' is larger"
		<< " than the size of the Index '" << index.getSize() << "'.",
		""
	);
	Index forward = index;
	forward[subindex]++;
	add(HoppingAmplitude(1/dx, index, index));
	add(HoppingAmplitude(-1/dx, index, forward));
}

};	//End of namespace MyTBTK

#endif
