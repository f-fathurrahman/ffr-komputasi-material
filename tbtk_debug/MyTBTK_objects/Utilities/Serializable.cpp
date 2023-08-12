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

/** @file Serializable.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/Serializable.h"
#include "MyTBTK/MyTBTKMacros.h"

#include "MyTBTK/json.hpp"

using namespace std;
//using namespace nlohmann;

namespace MyTBTK{

bool Serializable::validate(
	const string &serialization,
	const std::string &id,
	Mode mode
){
	switch(mode){
	case Mode::Debug:
	{
		size_t position = serialization.find(id + "(");
		if(position != 0 || serialization.back() != ')')
			return false;
/*		MyTBTKAssert(
			position == 0,
			"Serializable::validate()",
			"Unable to parse serialization string '"
			<< serialization << "'.",
			""
		);
		MyTBTKAssert(
			serialization.back() == ')',
			"Serializable::validate()",
			"Unable to parse serialization string '"
			<< serialization << "'.",
			""
		);*/

		return true;
	}
	case Mode::JSON:
	{
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			if(j.at("id").get<string>().compare(id) == 0)
				return true;
			else
				return false;
		}
		catch(nlohmann::json::exception &e){
			return false;
		}
	}
	default:
		MyTBTKExit(
			"Serializable::validate()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

bool Serializable::hasID(const string &serialization, Mode mode){
	switch(mode){
	case Mode::Debug:
	{
		size_t position = serialization.find('(');
		if(position == string::npos)
			return false;
		else
			return true;
	}
	case Mode::JSON:
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			j.at("id");
			return true;
		}
		catch(nlohmann::json::exception &e){
			return false;
		}
	default:
		MyTBTKExit(
			"Serializable::hasID()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

string Serializable::getID(const string &serialization, Mode mode){
	switch(mode){
	case Mode::Debug:
	{
		size_t position = serialization.find('(');
		MyTBTKAssert(
			position != string::npos,
			"Serializable::getID()",
			"Unable to parse serialization string '"
			<< serialization << "'.",
			""
		);

		return serialization.substr(0, position);
	}
	case Mode::JSON:
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			return j.at("id");
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"Serializable::getID()",
				"Unable to parse serialization string '"
				<< serialization << "'.",
				""
			);
		}
	default:
		MyTBTKExit(
			"Serializable::getID()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

string Serializable::extractComponent(
	const string &serialization,
	const string &containerID,
	const string &componentID,
	const string &componentName,
	Mode mode
){
	MyTBTKAssert(
		validate(serialization, containerID, mode),
		"Serializable::extractComponent()",
		"Invalid container ID in the serialization string '"
		<< serialization << "'. Expected '" << containerID << "'.",
		""
	);
	switch(mode){
	case Mode::JSON:
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);
			return j.at(componentName).dump();
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"Serializable::extractComponent()",
				"Unable to extract component with ID '"
				<< componentID << "' from the serialization"
				<< " string '" << serialization << "'.",
				""
			);
		}
	default:
		MyTBTKExit(
			"Serializable::extractComponent()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

string Serializable::getContent(const string &serialization, Mode mode){
	switch(mode){
	case Mode::Debug:
	{
		size_t position = serialization.find('(');
		MyTBTKAssert(
			position != string::npos,
			"Serializable::getContent()",
			"Unable to parse serialization string '"
			<< serialization << "'.",
			""
		);
		MyTBTKAssert(
			serialization.back() == ')',
			"Serializable::getContent()",
			"Unable to parse serialization string '"
			<< serialization << "'.",
			""
		);

		size_t contentStart = position + 1;
		size_t contentLength = serialization.size() - contentStart - 1;

		return serialization.substr(contentStart, contentLength);
	}
	default:
		MyTBTKExit(
			"Serializable::getContent()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

vector<string> Serializable::split(const string &content, Mode mode){
	switch(mode){
	case Mode::Debug:
	{
		vector<string> result;
		result.push_back(string());
		unsigned int bracketCounter = 0;
		for(unsigned int n = 0; n < content.size(); n++){
			char c = content.at(n);
			if(c == '(')
				bracketCounter++;
			if(c == ')'){
				MyTBTKAssert(
					bracketCounter > 0,
					"Serializable::split()",
					"Unable to split malformated serialization"
					<< " conent string '" << content << "'.",
					"Unbalanced brackets."
				);

				bracketCounter--;
			}

			if(c == ',' && bracketCounter == 0)
				result.push_back(string());
			else
				result.back() += c;
		}

		MyTBTKAssert(
			bracketCounter == 0,
			"Serializable::split()",
			"Unable to split malformated serialization conent string '"
			<< content << "'.",
			"Unbalanced brackets."
		);

		return result;
	}
	default:
		MyTBTKExit(
			"Serializable::split()",
			"Only Serializable::Mode::Debug is supported yet.",
			""
		);
	}
}

string Serializable::extract(
	const string &serialization,
	Mode mode,
	string component
){
	switch(mode){
	case Mode::JSON:
		try{
			nlohmann::json j = nlohmann::json::parse(serialization);

			return j.at(component).dump();
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"Serializable::extract()",
				"Unable to extract '" << component << "' from"
				<< " serialization string '" << serialization
				<< "'.",
				""
			);
		}
	default:
		MyTBTKExit(
			"Serializable::extract()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

};	//End of namespace MyTBTK
