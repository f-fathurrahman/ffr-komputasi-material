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

/** @file TransmissionRate.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/Property/TransmissionRate.h"
#include "MyTBTK/Streams.h"

#include "MyTBTK/json.hpp"

using namespace std;

namespace MyTBTK{
namespace Property{

TransmissionRate::TransmissionRate(
	double lowerBound,
	double upperBound,
	int resolution
) :
	EnergyResolvedProperty<double>(lowerBound, upperBound, resolution)
{
}

TransmissionRate::TransmissionRate(
	double lowerBound,
	double upperBound,
	int resolution,
	const double *data
) :
	EnergyResolvedProperty(lowerBound, upperBound, resolution, data)
{
}

TransmissionRate::TransmissionRate(
	const string &serialization, Mode mode
) :
	EnergyResolvedProperty(
		Serializable::extract(
			serialization,
			mode,
			"energyResolvedProperty"
		),
		mode
	)
{
	MyTBTKAssert(
		validate(serialization, "TransmissionRate", mode),
		"TransmissionRate::TransmissionRate()",
		"Unable to parse string as TransmissionRate '" << serialization
		<< "'.",
		""
	);

	switch(mode){
	case Mode::JSON:
		break;
	default:
		MyTBTKExit(
			"TransmissionRate::TransmissionRate()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

TransmissionRate& TransmissionRate::operator*=(const TransmissionRate &rhs){
	MyTBTKAssert(
		energyWindowsAreEqual(rhs),
		"TransmissionRate::operator*=()",
		"Incompatible energy windows.",
		""
	);

	for(unsigned int n = 0; n < getNumEnergies(); n++)
		operator()(n) *= rhs(n);

	return *this;
}

TransmissionRate& TransmissionRate::operator/=(const TransmissionRate &rhs){
	MyTBTKAssert(
		energyWindowsAreEqual(rhs),
		"TransmissionRate::operator/=()",
		"Incompatible energy windows.",
		""
	);

	for(unsigned int n = 0; n < getNumEnergies(); n++)
		operator()(n) /= rhs(n);

	return *this;
}

string TransmissionRate::serialize(Mode mode) const{
	switch(mode){
	case Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "TransmissionRate";
		j["energyResolvedProperty"] = nlohmann::json::parse(
			EnergyResolvedProperty::serialize(mode)
		);

		return j.dump();
	}
	default:
		MyTBTKExit(
			"TransmissionRate::serialize()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

};	//End of namespace Property
};	//End of namespace MyTBTK
