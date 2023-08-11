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

/// @cond MyTBTK_FULL_DOCUMENTATION
/** @package MyTBTKcalc
 *  @file SpacePartition.h
 *  @brief Base class for space partitions.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_SPACE_PARTITION
#define COM_DAFER45_MyTBTK_SPACE_PARTITION

#include "MyTBTK/Index.h"
#include "MyTBTK/Vector3d.h"

#include <vector>

namespace MyTBTK{

/** SpacePartition. Major indices indexes the unit cells, while minor indices
 *  arise when the unit cells are further subdivided. */
class SpacePartition{
public:
	/** Enum class used to indicate whether the mesh should be nodal or
	 *  interior. A nodal mesh is a mesh resulting from dividing a space
	 *  partition into equispaced elements and identifying the mesh with
	 *  the nodal grid points. An interior mesh is a mesh that results from
	 *  a similar division but where the mesh is identified with the
	 *  central point of the line/area/volume elements of the grid. Note
	 *  that interior points never fall on the bounding surface of the
	 *  space partition, while nodal points does. */
	enum class MeshType {Nodal, Interior};

	/** Constructor. */
	SpacePartition();

	/** Constructor. */
	SpacePartition(
		const std::vector<std::vector<double>> &basisVectors,
		MeshType meshType
	);

	/** Destructor. */
	virtual ~SpacePartition();

	/** Get number of dimensions. */
	unsigned int getNumDimensions() const;

	/** Returns the major index of the space partition for the given
	 *  coordinate. */
	virtual Index getMajorCellIndex(
		const std::vector<double> &coordinate
	) const = 0;

	/** Returns the minor index of the space partition corresponding to the
	 *  given coordinate, where the unit cell has been subdivided into
	 *  smaller cells specified by numMeshPoints. */
	virtual Index getMinorCellIndex(
		const std::vector<double> &coordinate,
		const std::vector<unsigned int> &numMeshPoints
	) const = 0;

	/** Returns an equispaced mesh spanned by the basis vectors, using
	 *  numMeshPoints mesh points along the corresponding directions. */
	virtual std::vector<std::vector<double>> getMajorMesh(
		const std::vector<unsigned int> &numMeshPoints
	) const = 0;

	/** Returns an equispaced mesh covering the unit cell, using
	 *  numMeshPoints mesh points along the corresponding directions. */
	virtual std::vector<std::vector<double>> getMinorMesh(
		const std::vector<unsigned int> &numMeshPoints
	) const = 0;

	/** Returns a single point of the minor mesh covering the unit cell,
	 *  using numMeshPoints mesh points along the corresponding directions.
	 *
	 *  @param meshPoint The mesh point to retrieve.
	 *  @param numMeshPoints. The number of mesh points in the minor mesh.
	 */
	virtual std::vector<double> getMinorMeshPoint(
		const std::vector<unsigned int> &meshPoint,
		const std::vector<unsigned int> &numMeshPoints
	) const = 0;
protected:
	/** Get basis vectors. */
	const std::vector<Vector3d>& getBasisVectors() const;

	/** Get mesh type. */
	MeshType getMeshType() const;
private:
	/** Lattice dimension. */
	unsigned dimensions;

	/** Basis vectors stored as three Vector3d. For lower dimensions the
	 *  vectors are padded with zeros to create three-dimensional vectors
	 *  and additional vectors along the extra dimensions are setup to unify
	 *  the calculations for all dimensions between 1-3. */
	std::vector<Vector3d> basisVectors;

	/** Mesh type. */
	MeshType meshType;
};

inline unsigned int SpacePartition::getNumDimensions() const{
	return dimensions;
}

inline const std::vector<Vector3d>& SpacePartition::getBasisVectors() const{
	return basisVectors;
}

inline SpacePartition::MeshType SpacePartition::getMeshType() const{
	return meshType;
}

};	//End namespace MyTBTK

#endif
/// @endcond
