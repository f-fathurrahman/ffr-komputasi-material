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

/** @package MyTBTKcalc
 *  @file Boolean.h
 *  @brief Boolean number.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_BOOLEAN
#define COM_DAFER45_MyTBTK_BOOLEAN

#include "MyTBTK/PseudoSerializable.h"

#ifndef MyTBTK_DISABLE_NLOHMANN_JSON
#	include "MyTBTK/json.hpp"
#endif

namespace MyTBTK{

#if MyTBTK_WRAP_PRIMITIVE_TYPES

/** @brief Boolean number. */
class Boolean : PseudoSerializable{
public:
	/** Constructor. */
	Boolean(){};

	/** Constructor.
	 *
	 *  @param value The value to initilize the Boolean with. */
	constexpr Boolean(bool value) : value(value) {}

	/** Constructs an Boolean from a serialization string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the Boolean.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	Boolean(const std::string &serialization, Serializable::Mode mode);

	/** Type conversion operator. */
	constexpr operator bool() const{	return value;	};

	/** Assignment operator.
	 *
	 *  @param value The value to assign the Boolean.
	 *
	 *  @return The Boolean after assignment has occured. */
	Boolean& operator=(bool rhs){
		value = rhs;

		return *this;
	}

	/** Serialize Boolean. Note that Boolean is PseudoSerializable rather
	 *  than Serializable. This means that the Serializable interface is
	 *  implemented non-virtually.
	 *
	 *  @param mode Serialization mode.
	 *
	 *  @return Serialized string representation of the Boolean. */
	std::string serialize(Serializable::Mode mode) const;

#ifndef MyTBTK_DISABLE_NLOHMANN_JSON
	/** Implements the Nlohmann json interface for conversion to json.
	 *
	 *  @param j The json output.
	 *  @param boolean The Boolean to convert. */
	friend void to_json(nlohmann::json &j, const Boolean &boolean){
		nlohmann::to_json(j, boolean.value);
	}

	/** Implements the Nlohmann json interface for conversion from json.
	 *
	 *  @param j The json input.
	 *  @param boolean The Boolean to convert to. */
	friend void from_json(const nlohmann::json &j, Boolean &boolean){
		nlohmann::from_json(j, boolean.value);
	}
#endif
private:
	/** Value. */
	bool value;
};

inline Boolean::Boolean(
	const std::string &serialization,
	Serializable::Mode mode
){
	switch(mode){
	case Serializable::Mode::JSON:
		value = stoi(serialization);
		break;
	default:
		MyTBTKExit(
			"Boolean::Boolean()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

inline std::string Boolean::serialize(Serializable::Mode mode) const{
	switch(mode){
	case Serializable::Mode::JSON:
		return std::to_string(value);
	default:
		MyTBTKExit(
			"Boolean::serialize()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

#else
	typedef int Boolean;
#endif

};	//End of namespace MyTBTK

#endif
