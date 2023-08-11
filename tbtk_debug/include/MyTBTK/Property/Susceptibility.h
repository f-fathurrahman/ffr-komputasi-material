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
 *  @file Susceptibility.h
 *  @brief Property container for the Susceptiility.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_PROPERTY_SUSCEPTIBILITY
#define COM_DAFER45_MyTBTK_PROPERTY_SUSCEPTIBILITY

#include "MyTBTK/Property/EnergyResolvedProperty.h"
#include "MyTBTK/MyTBTKMacros.h"

#include <complex>
#include <vector>

namespace MyTBTK{
namespace Property{

/** @brief Property container for the Susceptibility. */
class Susceptibility : public EnergyResolvedProperty<std::complex<double>>{
public:
	/** Constructs an uninitialized Susceptibility. */
	Susceptibility();

	/** Constructs a Susceptibility with real energies on the Custom
	 *  format. [See AbstractProperty for detailed information about the
	 *  Custom format.]
	 *
	 *  @param indexTree IndexTree containing the @link Index Indices
	 *  @endlink for which the Susceptibility should be contained.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy. */
	Susceptibility(
		const IndexTree &indexTree,
		double lowerBound,
		double upperBound,
		unsigned int resolution
	);

	/** Constructs a Susceptibility with real energies on the Custom format
	 *  and initializes it with data. [See AbstractProperty for detailed
	 *  information about the Custom format and the raw data format.]
	 *
	 *  @param indexTree IndexTree containing the @link Index Indices
	 *  @endlink for which the Susceptibility should be contained.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy.
	 *  @param data Raw data to initialize the GreensFunction with. */
	Susceptibility(
		const IndexTree &indexTree,
		double lowerBound,
		double upperBound,
		unsigned int resolution,
		const std::complex<double> *data
	);

	/** Constructs a Susceptibility with Matsubara energies on the Custom
	 *  format. [See AbstractProperty for detailed information about the
	 *  Custom format.]
	 *
	 *  @param indexTree IndexTree containing the @link Index Indices
	 *  @endlink for which the Susceptibility should be contained.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy. */
	Susceptibility(
		const IndexTree &indexTree,
		int lowerMatsubaraEnergyIndex,
		int upperMatsubaraEnergyIndex,
		double fundamentalMatsubaraEnergy
	);

	/** Constructs a Susceptibility with Matsubara energies on the Custom
	 *  format and initializes it with data. [See AbstractProperty for
	 *  detailed information about the Custom format and the raw data
	 *  format.]
	 *
	 *  @param indexTree IndexTree containing the @link Index Indices
	 *  @endlink for which the Susceptibility should be contained.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy.
	 *  @param data Raw data to initialize the GreensFunction with. */
	Susceptibility(
		const IndexTree &indexTree,
		int lowerMatsubaraEnergyIndex,
		int upperMatsubaraEnergyIndex,
		double fundamentalMatsubaraEnergy,
		const std::complex<double> *data
	);

	/** Constructor. Constructs the Susceptibility from a serialization
	 *  string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the Susceptibility.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	Susceptibility(const std::string &serialization, Mode mode);

	/** Overrides EnergyResolvedProperty::serialize(). */
	std::string serialize(Mode mode) const;
private:
};

};	//End namespace Property
};	//End namespace MyTBTK

#endif
/// @endcond
