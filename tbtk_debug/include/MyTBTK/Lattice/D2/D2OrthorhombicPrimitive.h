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
 *  @file D2OrthorhombicPrimitive.h
 *  @brief Orthorhombic primitive Bravais lattices.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_D2_ORTHORHOMBIC_PRIMITIVE
#define COM_DAFER45_MyTBTK_D2_ORTHORHOMBIC_PRIMITIVE

#include "MyTBTK/Lattice/D2/D2MonoclinicPrimitive.h"

namespace MyTBTK{
namespace Lattice{
namespace D2{

/** Orthorhombic primitiver Bravais lattice.
 *
 *  Dimensions:		2
 *  side0Length:	arbitrary
 *  side1Length:	arbitrary
 *  angle01:		pi/2 */
class OrthorhombicPrimitive : public MonoclinicPrimitive{
public:
	/** Constructor. */
	OrthorhombicPrimitive(
		double side0Length,
		double side1Length
	);

	/** Destructor. */
	~OrthorhombicPrimitive();
};

};	//End of namespace D2
};	//End of namespace Lattice
};	//End of namespace MyTBTK

#endif
/// @endcond
