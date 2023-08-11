/* Copyright 2016 Kristofer Björnson and Andreas Theiler
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
 *  @file FileParser.h
 *  @brief Reads and writes Model from and to text files.
 *
 *  @author Kristofer Björnson
 *  @author Andreas Theiler
 */

#ifndef COM_DAFER45_MyTBTK_FILE_PARSER
#define COM_DAFER45_MyTBTK_FILE_PARSER

#include "MyTBTK/Model.h"
#include "MyTBTK/ParameterSet.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace MyTBTK{

/** Reads and write
 */
class FileParser{
public:
	/** Enum for indicating storage mode. */
	enum class AmplitudeMode{
		ALL,
		ALL_EXCEPT_HC,
		UNIT_CELL,
		UNIT_CELL_EXCEPT_HC
	};

	/** Write Model to file. */
	static void writeModel(
		Model *model,
		std::string fileName,
		AmplitudeMode amplitudeMode,
		std::string description
	);

	/** Read Model from file. */
	static Model* readModel(std::string fileName);

	/** Write ParameterSet to file. */
	static void writeParameterSet(
		const ParameterSet &parameterSet,
		std::string fileName
	);

	/** Read ParameterSet from file. */
	static ParameterSet readParameterSet(std::string fileName);
private:
	/** Open output stream. */
	static void openOutput(std::string fileName);

	/** Close output stream. */
	static void closeOutput();

	/** Read input strem into internal input stream buffer (ssin). */
	static void readInput(std::string fileName);

	/** Write line breaks. */
	static void writeLineBreaks(int numLineBreaks);

	/** Write tabs. */
	static void writeTabs(int numTabs);

	/** Write a complex<double>. */
	static void write(std::complex<double> value);

	/** Write an Index. Example: [x y s]. */
	static void write(const Index &index);

	/** Write coordinates. Example: (0.1, 0.2, 0.3). */
	static void writeCoordinates(
//		const double *coordinates,
		const std::vector<double> &coordinates/*,
		int numCoordinates*/
	);

	/** Write specifiers. Example: <1 3>. */
//	static void writeSpecifiers(const int *specifiers, int numSpecifiers);
	static void writeSpecifiers(const std::vector<int> &specifiers);

	/** Write description comment. */
	static void writeDescription(std::string description);

	/** Write HoppingAmplitudes. */
	static void writeAmplitudes(Model *model, AmplitudeMode amplitudeMode);

	/** Write Geometry. */
	static void writeGeometry(Model *model);

	/** Reomve comments from file. */
	static void removeComments();

	/** Remove initial whitespaces. */
	static void removeInitialWhiteSpaces();

	/** Read a parameter */
	static int readParameter(
		std::string parameterName,
		std::string parentStructure
	);

	/** Read HoppingAmplitudes. */
	static void readAmplitudes(Model *model);

	/** Read Geometry. */
	static void readGeometry(Model *model);

	/** Read one HoppingAmplitude. */
	static HoppingAmplitude* readHoppingAmplitude();

	/** Read an Index. */
	static Index* readIndex();

	/** Read coordinates. Example (0.1, 0.2, 0.3). */
	static void readCoordinates(
		std::vector<double> *coordinates,
		int dimensions
	);

	/** Read specifiers. Example: <1 3>. */
	static void readSpecifiers(
		std::vector<int> *specifiers,
		int numSpecifiers
	);

	/** Read a complex<double>. */
	static bool readComplex(std::complex<double> *c);

	/** Read a double. */
	static bool readDouble(double *d, char endChar = ' ');

	/** Read an int. */
	static bool readInt(int *i, char endChar = ' ');

	/** Output file stream for writing. */
	static std::ofstream fout;

	/** String stream for parsing input. */
	static std::stringstream ssin;
};

};	//End of namespace MyTBTK

#endif
