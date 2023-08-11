/* Copyright 2019 Kristofer Björnson
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
 *  @file Polynomial.h
 *  @brief Class for storing polynomial expressions.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PADE_APPROXIMATOR_CONTINUOUS_FRACTIONS
#define COM_DAFER45_MyTBTK_PADE_APPROXIMATOR_CONTINUOUS_FRACTIONS

#include "MyTBTK/ArbitraryPrecision/Complex.h"
#include "MyTBTK/MyTBTKMacros.h"
#include "MyTBTK/Polynomial.h"

namespace MyTBTK{

class PadeApproximatorContinuousFractions{
public:
	/** Calculate the Padé approximation for a sampled function.
	 *
	 *  @param values The values of the function at the sample points.
	 *  @param The arguments at which the function has been sampled.
	 *
	 *  @return The continous fraction polynomial. */
	Polynomial<
		ArbitraryPrecision::Complex,
		ArbitraryPrecision::Complex,
		int
	> approximate(
		const std::vector<ArbitraryPrecision::Complex> &values,
		const std::vector<ArbitraryPrecision::Complex> &arguments
	);
};

}; //End of namespace MyTBTK

#endif
/// @endcond
