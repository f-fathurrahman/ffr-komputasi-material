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

/** @file CubicFaceCentered.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/Lattice/D3/CubicFaceCentered.h"
#include "MyTBTK/Vector3d.h"

using namespace std;

namespace MyTBTK{
namespace Lattice{
namespace D3{

CubicFaceCentered::CubicFaceCentered(double side0Length) :
	CubicPrimitive(
		side0Length
	)
{
	const vector<vector<double>> &latticeVectors = getLatticeVectors();

	Vector3d v0(latticeVectors.at(0));
	Vector3d v1(latticeVectors.at(1));
	Vector3d v2(latticeVectors.at(2));

	vector<vector<double>> additionalSites;
	additionalSites.push_back(((v0 + v1)/2.).getStdVector());
	additionalSites.push_back(((v0 + v2)/2.).getStdVector());
	additionalSites.push_back(((v1 + v2)/2.).getStdVector());

	setAdditionalSites(additionalSites);
}

CubicFaceCentered::~CubicFaceCentered(){
}

void CubicFaceCentered::makePrimitive(){
	const vector<vector<double>> &additionalSites = getAdditionalSites();

	vector<vector<double>> newLatticeVectors;
	newLatticeVectors.push_back(additionalSites.at(0));
	newLatticeVectors.push_back(additionalSites.at(1));
	newLatticeVectors.push_back(additionalSites.at(2));

	setLatticeVectors(newLatticeVectors);
	setAdditionalSites(vector<vector<double>>());
}

};	//End of namespace D3
};	//End of namespace Lattice
};	//End of namespace MyTBTK
