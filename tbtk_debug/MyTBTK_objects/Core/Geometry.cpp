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

/** @Geometry.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/Geometry.h"
#include "MyTBTK/Streams.h"
#include "MyTBTK/MyTBTKMacros.h"

#include "MyTBTK/json.hpp"

using namespace std;

namespace MyTBTK{

Geometry::Geometry(
){
	dimensions = -1;
}

Geometry::Geometry(
	const string &serialization,
	Mode mode
){
	MyTBTKAssert(
		validate(serialization, "Geometry", mode),
		"Geometry::Geometry()",
		"Unable to parse string as Geometry '" << serialization
		<< "'.",
		""
	);

	switch(mode){
	case Mode::JSON:
	{
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			dimensions = j.at("dimensions").get<int>();
			coordinates = IndexedDataTree<SerializableVector<double>>(
				j.at("coordinates").dump(),
				mode
			);
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"Geometry::Geometry()",
				"Unable to parse string as Geometry '"
				<< serialization << "'.",
				""
			);
		}

		break;
	}
	default:
		MyTBTKExit(
			"Geometry::Geometry()",
			"Only Serializable::Mode:Debug is supported yet.",
			""
		);
	}
}

Geometry::~Geometry(
){
}

void Geometry::translate(const vector<double> &translation){
	MyTBTKAssert(
		(int)translation.size() == dimensions,
		"Geometry::translate()",
		"Incompatible dimensions. The coordinates have '" << dimensions
		<< "' dimensions, but the translation has ' "
		<< translation.size() << "' dimensions.",
		""
	);

	for(
		IndexedDataTree<SerializableVector<double>>::Iterator iterator
			= coordinates.begin();
		iterator != coordinates.end();
		++iterator
	){
		for(int n = 0; n < dimensions; n++)
			(*iterator)[n] += translation[n];
	}
}

string Geometry::serialize(Mode mode) const{
	switch(mode){
	case Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "Geometry";
		j["dimensions"] = dimensions;
		j["coordinates"] = nlohmann::json::parse(
			coordinates.serialize(mode)
		);

		return j.dump();
	}
	default:
		MyTBTKExit(
			"Geometry::Geometry()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

};	//End of namespace MyTBTK
