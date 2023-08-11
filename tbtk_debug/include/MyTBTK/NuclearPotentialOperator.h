/* Copyright 2018 Kristofer Björnson
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
 *  @file NuclearPotentialOperator.h
 *  @brief Nuclear potential operator.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_NUCLEAR_POTENTIAL_OPERATOR
#define COM_DAFER45_MyTBTK_NUCLEAR_POTENTIAL_OPERATOR

#include "MyTBTK/AbstractOperator.h"
#include "MyTBTK/Atom.h"

namespace MyTBTK{

class NuclearPotentialOperator : public AbstractOperator{
public:
	/** Constructor. */
	NuclearPotentialOperator(
		const Atom &nucleus,
		const Vector3d &position
	);

	/** Get the nucleus.
	 *
	 *  @return An Atom specifying the type of the nucleus. */
	const Atom& getNucleus() const;

	/** Get position.
	 *
	 *  @return The position of the nucleus. */
	const Vector3d& getPosition() const;
private:
	/** An Atom specifying the type of the nucleus. */
	Atom nucleus;

	/** The position of the nucleus. */
	Vector3d position;
};

inline NuclearPotentialOperator::NuclearPotentialOperator(
	const Atom &nucleus,
	const Vector3d &position
) :
	AbstractOperator(OperatorID::NuclearPotential),
	nucleus(nucleus),
	position(position)
{
}

inline const Atom& NuclearPotentialOperator::getNucleus() const{
	return nucleus;
}

inline const Vector3d& NuclearPotentialOperator::getPosition() const{
	return position;
}

};	//End of namespace MyTBTK

#endif
/// @endcond
