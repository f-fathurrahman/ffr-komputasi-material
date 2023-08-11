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
 *  @file MyTBTKMacros.h
 *  @brief Precompiler macros
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_MACRO
#define COM_DAFER45_MyTBTK_MACRO

#include "MyTBTK/Streams.h"

#include <chrono>
#include <cstring>
#include <ctime>
#include <sstream>
#include <string>

#define MyTBTK_VERSION_STRING std::to_string(MyTBTK_VERSION_MAJOR) + "." \
	+ std::to_string(MyTBTK_VERSION_MINOR) + "." \
	+ std::to_string(MyTBTK_VERSION_PATCH)

#define MyTBTK_ABOUT_STRING \
	"MyTBTK\n" \
	"Version:\t" + MyTBTK_VERSION_STRING + "\n" \
	"Git hash:\t" MyTBTK_VERSION_GIT_HASH

inline std::string MyTBTK_GET_CURRENT_TIME_STRING(){
	std::chrono::time_point<std::chrono::system_clock> timePoint
		= std::chrono::system_clock::now();
	std::time_t now = std::chrono::system_clock::to_time_t(timePoint);

	return std::ctime(&now);
}

#define MyTBTK_RUNTIME_CONTEXT_STRING \
	MyTBTK_ABOUT_STRING + "\n" \
	+ "Date:\t" + MyTBTK_GET_CURRENT_TIME_STRING()

#ifdef MyTBTKOptimize
	#define MyTBTKAssert(expression, function, message, hint)	;
	#define MyTBTKExceptionAssert(expression, exception);
	#define MyTBTKExit(function, message, hint) exit(1);
#else
	#define MyTBTKAssert(expression, function, message, hint)	\
		if(!(expression)){	\
			MyTBTK::Streams::err << "Error in " << function << "\n";	\
			MyTBTK::Streams::err << "\t" << message << "\n";	\
			std::stringstream hintStream;	\
			hintStream << hint;	\
			if(std::strcmp(hintStream.str().c_str(), "") != 0)	\
				MyTBTK::Streams::err << "\tHint: " << hint << "\n";	\
			MyTBTK::Streams::err << "\tWhere: " << __FILE__ << ", " << __LINE__ << "\n";	\
			if(MyTBTK::Streams::logIsOpen())	\
				MyTBTK::Streams::closeLog();	\
			exit(1);	\
		}

	#define MyTBTKExceptionAssert(expression, exception)	\
		if(!(expression))	\
			throw exception;

	#define MyTBTKExit(function, message, hint)	\
		MyTBTK::Streams::err << "Error in " << function << "\n";	\
		MyTBTK::Streams::err << "\t" << message << "\n";	\
		std::stringstream hintStream;	\
		hintStream << hint;	\
		if(std::strcmp(hintStream.str().c_str(), "") != 0)	\
			MyTBTK::Streams::err << "\tHint: " << hint << "\n";	\
		MyTBTK::Streams::err << "\tWhere: " << __FILE__ << ", " << __LINE__ << "\n";	\
		if(MyTBTK::Streams::logIsOpen())	\
			MyTBTK::Streams::closeLog();	\
		exit(1);
#endif

#define MyTBTKNotYetImplemented(function)	\
	MyTBTK::Streams::err << "Error in " << function << "\n";	\
	MyTBTK::Streams::err << "\tNot yet implemented.\n";	\
	MyTBTK::Streams::err << "\tWhere: " << __FILE__ << ", " << __LINE__ << "\n";	\
	if(MyTBTK::Streams::logIsOpen())	\
		MyTBTK::Streams::closeLog();	\
	exit(1);

#define MyTBTKReadableCodeBlock(code) ;

#define MyTBTKWhere std::string(__FILE__) + ", " + std::to_string(__LINE__)

#endif
