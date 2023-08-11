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

/** @package MyTBTKcalc
 *  @file SpinPolarizedLDOS.h
 *  @brief Property container for spin-polarized local density of states
 *    (spin-polarized LDOS).
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_SPIN_POLARIZED_LDOS
#define COM_DAFER45_MyTBTK_SPIN_POLARIZED_LDOS

#include "MyTBTK/Property/EnergyResolvedProperty.h"
#include "MyTBTK/SpinMatrix.h"

#include <complex>

namespace MyTBTK{
namespace Property{

/** @brief Property container for spin-polarized local density of states
 *    (spin-polarized LDOS).
 *
 *  The SpinPolarizedLDOS is an EnergyResolvedProperty with DataType
 *  SpinMatrix.
 *
 *  # Example
 *  \snippet Property/SpinPolarizedLDOS.cpp SpinPolarizedLDOS
 *  ## Output
 *  \snippet output/Property/SpinPolarizedLDOS.txt SpinPolarizedLDOS
 *  \image html output/Property/SpinPolarizedLDOS/figures/PropertySpinPolarizedLDOSSpinPolarizedLDOS.png */
class SpinPolarizedLDOS : public EnergyResolvedProperty<SpinMatrix>{
public:
	/** Constructs a SpinPolarizedDOS on the Ranges format. [See
	 *  AbstractProperty for detailed information about the Ranges format.]
	 *
	 *  @param ranges The upper limits (exclusive) for the corresponding
	 *  dimensions.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy. */
	SpinPolarizedLDOS(
		const std::vector<int> &ranges,
		double lowerBound,
		double upperBound,
		int resolution
	);

	/** Constructs a SpinPolarizedDOS on the Ranges format and initializes
	 *  it with data. [See AbstractProperty for detailed information about
	 *  the Ranges format and the raw data format.]
	 *
	 *  @param ranges The upper limits (exclusive) for the corresponding
	 *  dimensions.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy.
	 *  @param data Raw data to initialize the SpinPolarizedLDOS with. */
	SpinPolarizedLDOS(
		const std::vector<int> &ranges,
		double lowerBound,
		double upperBound,
		int resolution,
		const SpinMatrix *data
	);

	/** Constructs a SpinPolarizedDOS on the Custom format. [See
	 *  AbstractProperty for detailed information about the Custom format.]
	 *
	 *  @param indexTree IndexTree containing the @link Index Indices
	 *  @endlink for which the SpinPolarizedLDOS is contained.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy. */
	SpinPolarizedLDOS(
		const IndexTree &indexTree,
		double lowerBound,
		double upperBound,
		int resolution
	);

	/** Constructs a SpinPolarizedDOS on the Custom format and initialize
	 *  it with data. [See AbstractProperty for detailed information about
	 *  the Custom format and the raw data format.]
	 *
	 *  @param indexTree IndexTree containing the @link Index Indices
	 *  @endlink for which the SpinPolarizedLDOS is contained.
	 *
	 *  @param lowerBound Lower bound for the energy.
	 *  @param upperBound Upper bound for the energy.
	 *  @param resolution Number of points to use for the energy.
	 *  @param data Raw data to initialize the SpinPolarizedLDOS with. */
	SpinPolarizedLDOS(
		const IndexTree &indexTree,
		double lowerBound,
		double upperBound,
		int resolution,
		const SpinMatrix *data
	);

	/** Constructor. Construct the SpinPolarizedLDOS from a serialization
	 *  string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the SpinPolarizedLDOS.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	SpinPolarizedLDOS(const std::string &serialization, Mode mode);

	/** Overrides Serializable::toString(). */
	virtual std::string toString() const;

	/** Overrides AbstractProperty::serialize(). */
	std::string serialize(Mode mode) const;
private:
};

};	//End namespace Property
};	//End namespace MyTBTK

#endif
