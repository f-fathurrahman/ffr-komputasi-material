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
 *  @file SourceAmplitude.h
 *  @brief Source amplitude for equations with a source term.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_SOURCE_AMPLITUDE
#define COM_DAFER45_MyTBTK_SOURCE_AMPLITUDE

#include "MyTBTK/Index.h"
#include "MyTBTK/Serializable.h"

#include <complex>
#include <initializer_list>
#include <vector>

namespace MyTBTK{

/** @brief Source amplitude for equations with a source term.
 *
 *  The Source Amplitude is used to encode information about the source term
 *  \f$S\f$ in an equation such as \f$i\hbar\frac{\partial\Psi}{\partial t} =
 *  H\Psi + S\f$. */
class SourceAmplitude : public Serializable{
public:
	/** Abstract base class for callbacks that allow for delayed
	 *  determination of the SourceAmplitude's value. */
	class AmplitudeCallback{
	public:
		/** Function responsible for returning the value of the
		 *  SourceAmplitude for the given Index.
		 *
		 *  @param index The Index to determine the value of the
		 *  SourceAmplitude for.
		 *
		 *  @return The value of the SourceAmplitude for the given
		 *  Index. */
		virtual std::complex<double> getSourceAmplitude(
			const Index &index
		) const = 0;
	};

	/** Constructs a SourceAmplitude. */
	SourceAmplitude();

	/** Constructs a SourceAmplitude from a value and an Index.
	 *
	 *  @param amplitude The amplitude value.
	 *  @param index The index for which the SourceAmplitude is defined. */
	SourceAmplitude(std::complex<double> amplitude, Index index);

	/** Constructor. Takes an AmplitudeCallback rather than a paramater
	 *  value. The AmplitudeCallback has to be defined such that it returns
	 *  a value for the given index when called at run time.
	 *
	 *  @param amplitudeCallback An AmplitudeCallback that is able to
	 *  return a value when passed an index.
	 *
	 *  @param index The Index for which the SourceAmplitude is defined. */
	SourceAmplitude(
		const AmplitudeCallback &amplitudeCallback,
		Index index
	);

	/** Constructor. Constructs the SourceAmplitude from a serialization
	 *  string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the SourceAmplitude.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	SourceAmplitude(
		const std::string &serializeation,
		Serializable::Mode mode
	);

	/** Get the amplitude value \f$S_{i}\f$.
	 *
	 *  @return The value of the amplitude. */
	std::complex<double> getAmplitude() const;

	/** Get index.
	 *
	 *  @return The Index. */
	const Index& getIndex() const;

	/** Get string representation of the SourceAmplitude.
	 *
	 *  @return A string representation of the SourceAmplitude. */
	std::string toString() const;

	/** Serialize SourceAmplitude. Note that SourceAmplitude is
	 *  pseudo-Serializable in that it implements the Serializable
	 *  interface, but does so non-virtually.
	 *
	 *  @param mode Serialization mode to use.
	 *
	 *  @return Serialized string representation of the SourceAmplitude.
	 */
	std::string serialize(Serializable::Mode mode) const;

	/** Get size in bytes.
	 *
	 *  @return Memory size required to store the SourceAmplitude. */
	unsigned int getSizeInBytes() const;
private:
	/** Amplitude \f$S_{i}\f$. Will be used if amplitudeCallback is NULL.
	 */
	std::complex<double> amplitude;

	/** AmplitudeCallback for runtime evaluation of amplitudes. Will be
	 *  called if not a nullptr. */
	const AmplitudeCallback *amplitudeCallback;

	/** Index at which the source is defined. */
	Index index;

};

inline SourceAmplitude::SourceAmplitude() : amplitudeCallback(nullptr){
}

inline std::complex<double> SourceAmplitude::getAmplitude() const{
	if(amplitudeCallback)
		return amplitudeCallback->getSourceAmplitude(index);
	else
		return amplitude;
}

inline const Index& SourceAmplitude::getIndex() const{
	return index;
}

inline std::string SourceAmplitude::toString() const{
	std::string str;
	str += "("
			+ std::to_string(real(amplitude))
			+ ", " + std::to_string(imag(amplitude))
		+ ")"
		+ ", " + index.toString();

	return str;
}

inline unsigned int SourceAmplitude::getSizeInBytes() const{
	return sizeof(SourceAmplitude)
		- sizeof(index)
		+ index.getSizeInBytes();
}

};	//End of namespace MyTBTK

#endif
