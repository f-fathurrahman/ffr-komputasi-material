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

/// @cond MyTBTK_FULL_DOCUMENTATION
/** @package MyTBTKcalc
 *  @file DefaultOperator.h
 *  @brief Default (dummy) operator class for indicating default behavior.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_DEFAULT_OPERATOR
#define COM_DAFER45_MyTBTK_DEFAULT_OPERATOR

#include "MyTBTK/AbstractOperator.h"

namespace MyTBTK{

class DefaultOperator : public AbstractOperator{
public:
	/** Constructor. */
	DefaultOperator();

	/** Destructor. */
	~DefaultOperator();
};

};	//End of namespace MyTBTK

#endif
/// @endcond
