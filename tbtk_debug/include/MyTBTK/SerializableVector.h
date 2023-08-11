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
 *  @file SerializableVector.h
 *  @brief Serializable wrapper of std::vector.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_SERIALIZABLE_VECTOR
#define COM_DAFER45_MyTBTK_SERIALIZABLE_VECTOR

#include "MyTBTK/Serializable.h"

#include <string>
#include <vector>

#ifndef MyTBTK_DISABLE_NLOHMANN_JSON
#	include "MyTBTK/json.hpp"
#endif

namespace MyTBTK{

template<typename DataType, bool = std::is_base_of<Serializable, DataType>::value>
class SerializableVector : public std::vector<DataType>, Serializable{
public:
	/** Constructor. */
	SerializableVector() : std::vector<DataType>(){};

	/** Constructor. */
	SerializableVector(
		const std::vector<DataType> &v
	) : std::vector<DataType>(v){};

	/** Constructor. Constructs the SerializableVector from a
	 *  serialization string.*/
	SerializableVector(const std::string &serialization, Mode mode);

	/** Implements Serializable::serialize(). */
	std::string serialize(Mode mode) const;
private:
};

template<typename DataType>
class SerializableVector<DataType, true> :
	public std::vector<DataType>, Serializable
{
public:
	/** Constructor. */
	SerializableVector() : std::vector<DataType>(){};

	/** Constructor. */
	SerializableVector(
		const std::vector<DataType> &v
	) : std::vector<DataType>(v){};

	/** Constructor. Constructs the SerializableVector from a
	 *  serialization string.*/
	SerializableVector(const std::string &serialization, Mode mode);

	/** Implements Serializable::serialize(). */
	std::string serialize(Mode mode) const;
private:
};

template<typename DataType>
class SerializableVector<DataType, false> :
	public std::vector<DataType>, Serializable
{
public:
	/** Constructor. */
	SerializableVector() : std::vector<DataType>(){};

	/** Constructor. */
	SerializableVector(
		const std::vector<DataType> &v
	) : std::vector<DataType>(v){};

	/** Constructor. Constructs the SerializableVector from a
	 *  serialization string.*/
	SerializableVector(const std::string &serialization, Mode mode);

	/** Implements Serializable::serialize(). */
	std::string serialize(Mode mode) const;
private:
};

#ifndef MyTBTK_DISABLE_NLOHMANN_JSON

template<typename DataType>
SerializableVector<DataType, false>::SerializableVector(
	const std::string &serialization,
	Mode mode
){
	MyTBTKAssert(
		validate(serialization, "SerializableVector", mode),
		"SerializableVector::SerializableVector()",
		"Unable to parse string as SerializableVector '" << serialization << "'.",
		""
	);

	switch(mode){
	case Mode::JSON:
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			nlohmann::json elements = j.at("elements");
			for(
				nlohmann::json::iterator it = elements.begin();
				it < elements.end();
				++it
			){
				std::vector<DataType>::push_back(
					Serializable::deserialize<DataType>(*it, mode)
				);
			}
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"SerializableVector::SerializableVector()",
				"Unable to parse string as SerializableVector"
				<< " '" << serialization << "'.",
				""
			);
		}

		break;
	default:
		MyTBTKExit(
			"SerializableVector::SerializableVector()",
			"Only SerializableVector::Mode::JSON is supported yet.",
			""
		);
	}
}

template<typename DataType>
SerializableVector<DataType, true>::SerializableVector(
	const std::string &serialization,
	Mode mode
){
	MyTBTKAssert(
		validate(serialization, "SerializableVector", mode),
		"SerializableVector::SerializableVector()",
		"Unable to parse string as SerializableVector '" << serialization << "'.",
		""
	);

	switch(mode){
	case Mode::JSON:
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			nlohmann::json elements = j.at("elements");
			for(
				nlohmann::json::iterator it = elements.begin();
				it < elements.end();
				++it
			){
				std::vector<DataType>::push_back(
					DataType(*it, mode)
				);
			}
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"SerializableVector::SerializableVector()",
				"Unable to parse string as SerializableVector"
				<< " '" << serialization << "'.",
				""
			);
		}

		break;
	default:
		MyTBTKExit(
			"SerializableVector::SerializableVector()",
			"Only SerializableVector::Mode::JSON is supported yet.",
			""
		);
	}
}

template<typename DataType>
std::string SerializableVector<DataType, false>::serialize(Mode mode) const{
	switch(mode){
	case Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "SerializableVector";
		j["elements"] = nlohmann::json::array();
		for(unsigned int n = 0; n < std::vector<DataType>::size(); n++){
			j["elements"].push_back(
				Serializable::serialize(
					std::vector<DataType>::at(n),
					mode
				)
			);
		}

		return j.dump();
	}
	default:
		MyTBTKExit(
			"SerializableVector::serialize()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

template<typename DataType>
std::string SerializableVector<DataType, true>::serialize(Mode mode) const{
	switch(mode){
	case Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "SerializableVector";
		j["elements"] = nlohmann::json::array();
		for(unsigned int n = 0; n < std::vector<DataType>::size(); n++){
			j["elements"].push_back(
				std::vector<DataType>::at(n).serialize(mode)
			);
		}

		return j.dump();
	}
	default:
		MyTBTKExit(
			"SerializableVector::serialize()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

#endif

};	//End of namespace MyTBTK

#endif
/// @endcond
