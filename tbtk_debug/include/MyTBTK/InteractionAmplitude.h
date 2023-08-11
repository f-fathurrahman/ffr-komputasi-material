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
 *  @file InteractionAmplitude.h
 *  @brief Interaction amplitude.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_INTERACTION_AMPLITUDE
#define COM_DAFER45_MyTBTK_INTERACTION_AMPLITUDE

#include "MyTBTK/Index.h"

#include <complex>
#include <initializer_list>

namespace MyTBTK{

class InteractionAmplitude{
public:
	/** Constructor. */
	InteractionAmplitude();

	/** Constructor. */
	InteractionAmplitude(
		std::complex<double> amplitude,
		std::initializer_list<Index> creationOperatorIndices,
		std::initializer_list<Index> annihilationOperatorIndices
	);

	/** Constructor. Takes a callback function rather than a parameter
	 *  value. The callback function has to be defined such that it returns
	 *  a value for the given indices when called at run time. */
	InteractionAmplitude(
		std::complex<double> (*amplitudeCallback)(const std::vector<Index>&, const std::vector<Index>&),
		std::initializer_list<Index> toIndex,
		std::initializer_list<Index> fromIndex
	);

	/** Copy constructor. */
	InteractionAmplitude(const InteractionAmplitude &ia);

	/** Destructor. */
	~InteractionAmplitude();

	/** Get the amplitude value \f$a_{\{c\}\{a\})}\f$, where \f$\{c\}\f$
	 *  and \f$\{a\}\f$ are sets of creation and annihilation operator
	 *  indices. */
	std::complex<double> getAmplitude() const;

	/** Returns the number of creation operator indices. */
	unsigned int getNumCreationOperators() const;

	/** Returns the number of annihilation operator indices. */
	unsigned int getNumAnnihilationOperators() const;

	/** Get creation operator index. */
	const Index& getCreationOperatorIndex(unsigned int n) const;

	/** Get annihilation operator index. */
	const Index& getAnnihilationOperatorIndex(unsigned int n) const;
private:
	/** Amplitude. */
	std::complex<double> amplitude;

	/** Callback function for runtime evaluation of amplitudes. Will be
	 *  called if not null. */
	std::complex<double> (*amplitudeCallback)(
		const std::vector<Index> &toIndex,
		const std::vector<Index> &fromIndex
	);

	/** Indices for creation operators. */
	std::vector<Index> creationOperatorIndices;

	/** Indices for annihilation operators. */
	std::vector<Index> annihilationOperatorIndices;
};

inline std::complex<double> InteractionAmplitude::getAmplitude() const{
	if(amplitudeCallback)
		return amplitudeCallback(creationOperatorIndices, annihilationOperatorIndices);
	else
		return amplitude;
}

inline unsigned int InteractionAmplitude::getNumCreationOperators() const{
	return creationOperatorIndices.size();
}

inline unsigned int InteractionAmplitude::getNumAnnihilationOperators() const{
	return annihilationOperatorIndices.size();
}

inline const Index& InteractionAmplitude::getCreationOperatorIndex(unsigned int n) const{
	return creationOperatorIndices.at(n);
}

inline const Index& InteractionAmplitude::getAnnihilationOperatorIndex(unsigned int n) const{
	return annihilationOperatorIndices.at(n);
}

};	//End of namespace MyTBTK

#endif
/// @endcond
