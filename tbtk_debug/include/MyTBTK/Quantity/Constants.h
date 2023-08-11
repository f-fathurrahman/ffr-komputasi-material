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

/// @cond MyTBTK_FULL_DOCUMENTATION
/** @package MyTBTKcalc
 *  @file Constants.h
 *  @brief Numerical constants.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_QUANTITY_CONSTANTS
#define COM_DAFER45_MyTBTK_QUANTITY_CONSTANTS

#include "MyTBTK/Quantity/Base.h"
#include "MyTBTK/Quantity/Constant.h"
#include "MyTBTK/Quantity/Derived.h"
#include "MyTBTK/MyTBTK.h"
#include "MyTBTK/UnitHandler.h"

#include <cmath>
#include <map>
#include <string>

namespace MyTBTK{
namespace Quantity{

/** @brief Numerical constants.
 *
 *  <b>Warning:</b> All numerical constants defined here are given in the
 *  default base units C, pcs, eV, m, K, s. These are used to initialize the
 *  UnitHandler and Quantities and should not be used directly outside of this
 *  scope. Therefore, request constants through the UnitHandler for all other
 *  purposes. */
class Constants{
public:
	/** Get a constant with a given name.
	 *
	 *  @param name The name of the constant. */
	static Constant get(const std::string &name);
private:
	/** Container for the constants. */
	static std::map<std::string, Constant> constants;

	/** Function to be called from MyTBTKs global initialization function. */
	static void initialize();

	/** The global initialization function is a friend to allow for
	 *  initialization of the constants. */
	friend void MyTBTK::Initialize();
};

inline Constant Constants::get(const std::string &name){
	try{
		return constants.at(name);
	}
	catch(const std::out_of_range &e){
		MyTBTKExit(
			"Constants::get()",
			"Unknown constant '" << name << "'.",
			""
		);
	}
}

}; //End of namesapce Quantity
}; //End of namesapce MyTBTK

#endif
/// @endcond
