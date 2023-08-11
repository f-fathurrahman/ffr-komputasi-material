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

/** @file SingleParticleContext.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/SingleParticleContext.h"

#include "MyTBTK/json.hpp"

using namespace std;

namespace MyTBTK{

SingleParticleContext::SingleParticleContext(){
	statistics = Statistics::FermiDirac;
}

SingleParticleContext::SingleParticleContext(
	const vector<unsigned int> &capacity
) :
	hoppingAmplitudeSet(capacity)
{
	statistics = Statistics::FermiDirac;
}

SingleParticleContext::SingleParticleContext(
	const string &serialization,
	Mode mode
)
{
	MyTBTKAssert(
		validate(serialization, "SingleParticleContext", mode),
		"SingleParticleContext::SingleParticleContext()",
		"Unable to parse string as SingleParticleContext '"
		<< serialization << "'.",
		""
	);

	switch(mode){
	case Mode::JSON:
	{
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			statistics = Serializable::deserialize<Statistics>(
				j.at("statistics").get<string>(),
				mode
			);
			hoppingAmplitudeSet = HoppingAmplitudeSet(
				j.at("hoppingAmplitudeSet").dump(),
				mode
			);
			geometry = Geometry(
				j.at("geometry").dump(),
				mode
			);
			sourceAmplitudeSet = SourceAmplitudeSet(
				j.at("sourceAmplitudeSet").dump(),
				mode
			);
			overlapAmplitudeSet = OverlapAmplitudeSet(
				j.at("overlapAmplitudeSet").dump(),
				mode
			);
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"SingleParticleContext::SingleParticleContext()",
				"Unable to parse string as"
				<< " SingleParticleContext '" << serialization
				<< "'.",
				""
			);
		}

		break;
	}
	default:
		MyTBTKExit(
			"SingleParticleContext::SingleParticleContext()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

void SingleParticleContext::generateHoppingAmplitudeSet(
	const HoppingAmplitude::AmplitudeCallback &hoppingAmplitudeCallback
){
	for(
		BasisStateSet::ConstIterator toIterator
			= basisStateSet.cbegin();
		toIterator != basisStateSet.cend();
		++toIterator
	){
		for(
			BasisStateSet::ConstIterator fromIterator
				= basisStateSet.cbegin();
			fromIterator != basisStateSet.cend();
			++fromIterator
		){
			hoppingAmplitudeSet.add(
				HoppingAmplitude(
					hoppingAmplitudeCallback,
					(*toIterator).getIndex(),
					(*fromIterator).getIndex()
				)
			);
		}
	}
}

void SingleParticleContext::generateOverlapAmplitudeSet(
	const OverlapAmplitude::AmplitudeCallback &overlapAmplitudeCallback
){
	for(
		BasisStateSet::ConstIterator toIterator
			= basisStateSet.cbegin();
		toIterator != basisStateSet.cend();
		++toIterator
	){
		for(
			BasisStateSet::ConstIterator fromIterator
				= basisStateSet.cbegin();
			fromIterator != basisStateSet.cend();
			++fromIterator
		){
			overlapAmplitudeSet.add(
				OverlapAmplitude(
					overlapAmplitudeCallback,
					(*toIterator).getIndex(),
					(*fromIterator).getIndex()
				)
			);
		}
	}
}

string SingleParticleContext::serialize(Mode mode) const{
	switch(mode){
	case Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "SingleParticleContext";
		j["statistics"] = Serializable::serialize(statistics, mode);
		j["hoppingAmplitudeSet"] = nlohmann::json::parse(
			hoppingAmplitudeSet.serialize(mode)
		);
		j["geometry"] = nlohmann::json::parse(
			geometry.serialize(mode)
		);
		j["sourceAmplitudeSet"] = nlohmann::json::parse(
			sourceAmplitudeSet.serialize(mode)
		);
		j["overlapAmplitudeSet"] = nlohmann::json::parse(
			overlapAmplitudeSet.serialize(mode)
		);

		return j.dump();
	}
	default:
		MyTBTKExit(
			"SingleParticleContext::serialize()",
			"Only Serializable::Mode::Debugis supported yet.",
			""
		);
	}
}

};	//End of namespace MyTBTK
