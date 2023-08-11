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

/** @file HoppingAmplitude.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/HoppingAmplitude.h"
#include "MyTBTK/Streams.h"

#include <sstream>

#include "MyTBTK/json.hpp"

using namespace std;

namespace MyTBTK{

HoppingAmplitude::HoppingAmplitude(){
	this->amplitudeCallback = nullptr;
}

HoppingAmplitude::HoppingAmplitude(
	complex<double> amplitude,
	Index toIndex,
	Index fromIndex
) :
	fromIndex(fromIndex),
	toIndex(toIndex)
{
	this->amplitude = amplitude;
	this->amplitudeCallback = NULL;
};

HoppingAmplitude::HoppingAmplitude(
	const AmplitudeCallback &amplitudeCallback,
	Index toIndex,
	Index fromIndex
) :
	fromIndex(fromIndex),
	toIndex(toIndex)
{
	this->amplitudeCallback = &amplitudeCallback;
};

HoppingAmplitude::HoppingAmplitude(
	const HoppingAmplitude &ha
) :
	fromIndex(ha.fromIndex),
	toIndex(ha.toIndex)
{
	amplitude = ha.amplitude;
	this->amplitudeCallback = ha.amplitudeCallback;
}

HoppingAmplitude::HoppingAmplitude(
	const string &serialization,
	Serializable::Mode mode
){
	MyTBTKAssert(
		Serializable::validate(
			serialization,
			"HoppingAmplitude",
			mode
		),
		"HoppingAmplitude::HoppingAmplitude()",
		"Unable to parse string as HoppingAmplitude '"
		<< serialization << "'.",
		""
	);

	switch(mode){
	case Serializable::Mode::Debug:
	{
		string content = Serializable::getContent(serialization, mode);

		vector<string> elements = Serializable::split(
			content,
			Serializable::Mode::Debug
		);

		amplitudeCallback = nullptr;

		stringstream ss;
		ss.str(elements.at(0));
		ss >> amplitude;
		toIndex = Index(elements.at(1), mode);
		fromIndex = Index(elements.at(2), mode);
		break;
	}
	case Serializable::Mode::JSON:
	{
		try{
			amplitudeCallback = nullptr;

			nlohmann::json j = nlohmann::json::parse(serialization);
//			Serializable::deserialize(j["amplitude"].get<string>(), &amplitude, mode);
			amplitude = Serializable::deserialize<complex<double>>(
				j["amplitude"].get<string>(),
				mode
			);
			toIndex = Index(j["toIndex"].dump(), mode);
			fromIndex = Index(j["fromIndex"].dump(), mode);
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"HoppingAmplitude::HoppingAmplitude()",
				"Unable to parse string as HoppingAmplitude '"
				<< serialization << "'.",
				""
			);
		}

		break;
	}
	default:
		MyTBTKExit(
			"HoppingAmplitude::HoppingAmplitude()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

HoppingAmplitude HoppingAmplitude::getHermitianConjugate() const{
	if(amplitudeCallback)
		return HoppingAmplitude(*amplitudeCallback, fromIndex, toIndex);
	else
		return HoppingAmplitude(conj(amplitude), fromIndex, toIndex);
}

void HoppingAmplitude::print() const{
	Streams::out << "From index:\t";
	for(unsigned int n = 0; n < fromIndex.getSize(); n++){
		Streams::out << fromIndex.at(n) << " ";
	}
	Streams::out << "\n";
	Streams::out << "To index:\t";
	for(unsigned int n = 0; n < toIndex.getSize(); n++){
		Streams::out << toIndex.at(n) << " ";
	}
	Streams::out << "\n";
	Streams::out << "Amplitude:\t" << getAmplitude() << "\n";
}

string HoppingAmplitude::serialize(Serializable::Mode mode) const{
	MyTBTKAssert(
		amplitudeCallback == nullptr,
		"HoppingAmplitude::serialize()",
		"Unable to serialize HoppingAmplitude that uses callback"
		<< " value.",
		""
	);

	switch(mode){
	case Serializable::Mode::Debug:
	{
		stringstream ss;
		ss << "HoppingAmplitude(";
		ss << Serializable::serialize(amplitude, mode) << ",";
		ss << toIndex.serialize(mode) << "," << fromIndex.serialize(mode);
		ss << ")";

		return ss.str();
	}
	case Serializable::Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "HoppingAmplitude";
		j["amplitude"] = Serializable::serialize(amplitude, mode);
		j["toIndex"] = nlohmann::json::parse(toIndex.serialize(mode));
		j["fromIndex"] = nlohmann::json::parse(fromIndex.serialize(mode));

		return j.dump();
	}
	default:
		MyTBTKExit(
			"HoppingAmplitude::serialize()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

};	//End of namespace MyTBTK
