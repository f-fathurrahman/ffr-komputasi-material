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
 *  @file Index.h
 *  @brief Physical index.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_INDEX
#define COM_DAFER45_MyTBTK_INDEX

#include "MyTBTK/Subindex.h"
#include "MyTBTK/Serializable.h"
#include "MyTBTK/Streams.h"

#include <vector>

namespace MyTBTK{

/** @brief Physical index.
 *
 *  An Index is a list of @link Subindex Subindices@endlink. For example,
 *  {x, y, spin}, {x, y, z, orbital, spin}, and
 *  {subsystem, x, y, z, orbital, spin}.
 *
 *  # Example
 *  \snippet Core/Index.cpp Index
 *  ## Output
 *  \snippet output/Core/Index.txt Index */
class Index{
public:
	//MyTBTKFeature Core.Index.Construction.1 2019-09-19
	/** Constructs an empty Index. */
	Index(){};

	//MyTBTKFeature Core.Index.Construction.2 2019-09-19
	/** Constructs an Index from an initializer list.
	 *
	 * @param i Initializer list from which the Index is constructed. */
	Index(std::initializer_list<Subindex> i) : indices(i){};

	/** Constructs an Index from an initializer list.
	 *
	 * @param i Initializer list from which the Index is constructed. */
/*	Index(
		std::initializer_list<unsigned int> i
	) : indices(i.begin(), i.end()){};*/

	//MyTBTKFeature Core.Index.Construction.3.C++ 2019-09-19
	/** Constructs an Index from an std::vector<int>.
	 *
	 *  @param i Vector from which the Index is constructed. */
	Index(std::vector<Subindex> i) : indices(i){};

	/** Constructs an Index from an std::vector<int>.
	 *
	 *  @param i Vector from which the Index is constructed. */
//	Index(std::vector<unsigned int> i) : indices(i.begin(), i.end()){};

	//MyTBTKFeature Core.Index.Copy.1 2019-09-19
	//MyTBTKFeature Core.Index.Copy.2 2019-09-23
	/** Copy constructor.
	 *
	 *  @param index Index to copy. */
	Index(const Index &index) : indices(index.indices){};

	//MyTBTKFeature Core.Index.Construction.4 2019-09-19
	//MyTBTKFeature Core.Index.Construction.5.C++ 2019-09-19
	/** Constructs a new Index by concatenating two indices into one total
	 *  index of the form {head, tail}.
	 *
	 *  @param head First part of the compund Index.
	 *  @param tail Second part of the compund Index.*/
	Index(const Index &head, const Index &tail);

	/** Constructs a compound Index by concatenating a list of indices,
	 *  adding IDX_SEPARATOR between every index.
	 *
	 *  @param indexList List of indices. */
//	Index(std::initializer_list<std::initializer_list<Subindex>> indexList);

	//MyTBTKFeature Core.Index.Construction.8.C++ 2019-09-19
	/** Constructs a compund Index by concatenating a list of indices,
	 *  adding IDX_SEPARATOR between every index.
	 *
	 *  @param indexList List of indices. */
	Index(const std::vector<std::vector<Subindex>> &indexList);

	//MyTBTKFeature Core.Index.Construction.6 2019-09-19
	//MyTBTKFeature Core.Index.Construction.7.C++ 2019-09-19
	/** Constructs a compound Index by Concatenating a list of indices,
	 *  adding IDX_SEPARATOR between every index.
	 *
	 *  @param indexList List of indices. */
	Index(std::initializer_list<Index> indexList);

	/** Constructs a compound Index by Concatenating a list of indices,
	 *  adding IDX_SEPARATOR between every index.
	 *
	 *  @param indexList List of indices. */
	Index(std::vector<Index> indexList);

	//MyTBTKFeature Core.Index.Construction.9 2019-09-19
	/** Constructs an Index from a string.
	 *
	 *  @param indexString String such as "{1, 2, 3} from which the Index is
	 *  constructed. */
	Index(const std::string &indexString);

	//MyTBTKFeature Core.Index.Serialization.1 2019-09-19
	/** Constructs an Index from a serialization string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the Index.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	Index(const std::string &serialization, Serializable::Mode mode);

	//MyTBTKFeature Core.Index.Equals.1 2019-09-19
	//MyTBTKFeature Core.Index.Equals.2 2019-09-19
	//MyTBTKFeature Core.Index.Equals.3 2019-09-19
	//MyTBTKFeature Core.Index.Equals.4 2019-09-21
	//MyTBTKFeature Core.Index.Equals.5 2019-09-19
	//MyTBTKFeature Core.Index.Equals.6 2019-09-19
	//MyTBTKFeature Core.Index.Equals.7 2019-09-19
	//MyTBTKFeature Core.Index.Equals.8 2019-09-19
	//MyTBTKFeature Core.Index.Equals.9 2019-09-19
	//MyTBTKFeature Core.Index.Equals.10 2019-09-19
	/** Compare this index with another index. Returns true if the indices
	 *  have the same number of subindices and all subindices are equal.
	 *
	 *  @param index Index to compare with.
	 *  @param allowWildcard IDX_ALL is interpreted as wildcard.
	 *
	 *  @return True if the indices are equal, otherwise false. */
	bool equals(const Index &index, bool allowWildcard = false) const;

	//MyTBTKFeature Core.Index.at.1 2019-09-19
	/** Get subindex n.
	 *
	 *  @param n Subindex.
	 *
	 *  @return Subindex at position n. */
	Subindex& at(unsigned int n);

	//MyTBTKFeature Core.Index.at.2.C++ 2019-09-19
	/** Get subindex n. Constant version.
	 *
	 *  @param n Subindex.
	 *
	 *  @return Subindex at position n. */
	const Subindex& at(unsigned int n) const;

	//MyTBTKFeature Core.Index.getSize.1 2019-09-19
	//MyTBTKFeature Core.Index.getSize.2 2019-09-19
	//MyTBTKFeature Core.Index.getSize.3 2019-09-19
	//MyTBTKFeature Core.Index.getSize.4 2019-09-19
	//MyTBTKFeature Core.Index.getSize.5 2019-09-19
	/** Get size.
	 *
	 *  @return Number of subindices for individual indices such as
	 *  {1, 2, 3}. For compound indices such as {{1, 2}, {3,4}, {5, 6}},
	 *  the total number of subindices (here 6) plus the number of index
	 *  separators (here 2) are returned (here 6+2=8). */
	unsigned int getSize() const;

	/** Reserves memory for the Index.
	 *
	 *  @param Number of subindices to reserve space for. */
	void reserve(unsigned int size);

	//MyTBTKFeature Core.Index.pushBack.1 2019-09-19
	/** Push subindex at the back of the index.
	 *
	 *  @param subindex Subindex to append to the Index. */
	void pushBack(Subindex subindex);

	//MyTBTKFeature Core.Index.popFront.1 2019-09-19
	/** Removes and returns the first subindex.
	 *
	 *  @return The first subindex. */
	Subindex popFront();

	//MyTBTKFeature Core.Index.popBack.1 2019-09-19
	/** Removes and returns the last subindex.
	 *
	 *  @return The last subindex. */
	Subindex popBack();

	//MyTBTKFeature Core.Index.insert.1 2019-10-28
	/** Insert a Subindex at a given position.
	 *
	 *  @param n Subindex position to insert at.
	 *  @param subindex Subindex to insert. */
	void insert(unsigned int n, Subindex subindex);

	//MyTBTKFeature Core.Index.erase.1 2019-10-28
	/** Remove and return the Subindex at a given position.
	 *
	 *  @param n Subindex position to remove.
	 *
	 *  @return The value of the removed Subindex. */
	Subindex erase(unsigned int n);

	//MyTBTKFeature Core.Index.getUnitRange.1 2019-09-19
	/** Returns an index with the same number or subindices, and each
	 *  subindex set to 1.
	 *
	 *  @return Index with all subindices set to 1. */
	Index getUnitRange();

	//MyTBTKFeature Core.Index.getSubIndex.1 2019-09-21
	/** Returns an Index containing the subindices from position 'first' to
	 *  'last'.
	 *
	 *  @parameter first First index to include in range (inclusive).
	 *  @parameter last Last index to include in range (inclusive).
	 *
	 *  @return An index containing the subindices in the range first to
	 *  last (inclusive). */
	Index getSubIndex(int first, int last) const;

	//MyTBTKFeature Core.Index.split.1 2019-09-19
	/** Split a compound Index into its components.
	 *
	 *  @return An std::vector<Index> containing the individual @link Index
	 *  Indices @endlink.*/
	std::vector<Index> split() const;

	//MyTBTKFeature Core.Index.isPatternIndex.1 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.2 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.3 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.4 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.5 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.6 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.7 2019-09-19
	//MyTBTKFeature Core.Index.isPatternIndex.8 2019-09-19
	/** Returns true if the Index is a pattern index. That is, if it
	 *  contains a negative subindex.
	 *
	 *  @return True if the Index is a pattern index, otherwise false. */
	bool isPatternIndex() const;

	/** Print index. Mainly for debuging. */
	void print() const;

	//MyTBTKFeature Core.Index.toString.1 2019-09-19
	//MyTBTKFeature Core.Index.toString.2 2019-09-19
	/** Get string representation of the Index.
	 *
	 *  @return A string representation of the Index. */
	std::string toString() const;

	/** Writes the Index toString()-representation to a stream.
	 *
	 *  @param stream The stream to write to.
	 *  @param index The Index to write.
	 *
	 *  @return Reference to the output stream just written to. */
	friend std::ostream& operator<<(
		std::ostream &stream,
		const Index &index
	);

	//MyTBTKFeature Core.Index.operator<.1.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator<.2.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator<.3.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator<.4.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator<.5.C++ 2019-09-19
	/** Comparison operator. Returns false if the TreeNode structure would
	 *  generate a smaller Hilbert space index for i1 than for i2.
	 *
	 *  @return True if i1 would generate a smaller Hilbert space index
	 *  than i2. */
	friend bool operator<(const Index &i1, const Index &i2);

	//MyTBTKFeature Core.Index.operator>.1.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator>.2.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator>.3.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator>.4.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator>.5.C++ 2019-09-19
	/** Comparison operator. Returns false if the TreeNode structure would
	 *  generate a larger Hilbert space index for i1 than for i2.
	 *
	 *  @return True if i1 would generate a larger Hilbert space index than
	 *  i2. */
	friend bool operator>(const Index &i1, const Index &i2);

	//MyTBTKFeature Core.Index.operator[].1.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator[].2.C++ 2019-09-19
	/** Subscript operator.
	 *
	 *  @param n Subindex.
	 *
	 *  @return Subindex at position n. */
	Subindex& operator[](unsigned int subindex);

	//MyTBTKFeature Core.Index.operator[].1.C++ 2019-09-19
	//MyTBTKFeature Core.Index.operator[].2.C++ 2019-09-19
	/** Subscript operator.
	 *
	 *  @param n Subindex.
	 *
	 *  @return Subindex at position n. */
	const Subindex& operator[](unsigned int subindex) const;

	//MyTBTKFeature Core.Index.Serialization.1 2019-09-19
	/** Serialize Index. Note that Index is pseudo-Serializable in that it
	 *  implements the Serializable interface, but does so non-virtually.
	 *
	 *  @param mode Serialization mode to use.
	 *
	 *  @return Serialized string represenation of the Index. */
	std::string serialize(Serializable::Mode mode) const;

	/** Get size in bytes.
	 *
	 *  @return Memory size required to store the Index. */
	unsigned int getSizeInBytes() const;
private:
	/** Subindex container. */
	std::vector<Subindex> indices;
};

inline void Index::print() const{
	Streams::out << "{";
	for(unsigned int n = 0; n < indices.size(); n++){
		if(n != 0)
			Streams::out << ", ";
		Streams::out << indices.at(n);
	}
	Streams::out << "}\n";
}

inline std::string Index::toString() const{
	std::string str = "{";
	bool isFirstIndex = true;
	for(unsigned int n = 0; n < indices.size(); n++){
		Subindex subindex = indices.at(n);
		if(!isFirstIndex && !subindex.isIndexSeparator())
			str += ", ";
		else
			isFirstIndex = false;
		if(subindex.isWildcard()){
			str += "IDX_ALL";
		}
		else if(subindex.isSummationIndex()){
			str += "IDX_SUM_ALL";
		}
		else if(subindex.isRangeIndex()){
			switch(subindex){
			case IDX_X:
				str += "IDX_X";
				break;
			case IDX_Y:
				str += "IDX_Y";
				break;
			case IDX_Z:
				str += "IDX_Z";
				break;
			default:
				MyTBTKExit(
					"Index::toString()",
					"This should never happen, contact the"
					<< " developer.",
					""
				);
			}
		}
		else if(subindex.isSpinIndex()){
			str += "IDX_SPIN";
		}
		else if(subindex.isIndexSeparator()){
			str += "}, {";
			isFirstIndex = true;
		}
		else if(subindex.isLabeledWildcard()){
			str += "IDX_ALL_(";
			str += std::to_string(subindex.getWildcardLabel());
			str += ")";
		}
		else{
			str += std::to_string(subindex);
		}
	}
	str += "}";

	return str;
}

inline std::ostream& operator<<(std::ostream &stream, const Index &index){
	stream << index.toString();

	return stream;
}

inline bool Index::equals(const Index &index, bool allowWildcard) const{
	if(indices.size() == index.indices.size()){
		for(unsigned int n = 0; n < indices.size(); n++){
			if(indices.at(n) != index.indices.at(n)){
				if(!allowWildcard)
					return false;
				else{
					if(
						indices.at(n).isWildcard() ||
						index.indices.at(n).isWildcard()
					){
						continue;
					}
					else if(
						indices.at(n).isLabeledWildcard()
					){
						for(
							unsigned int c = 0;
							c < indices.size();
							c++
						){
							if(
								indices.at(c)
								== indices.at(n)
								&& index.indices.at(c)
								!= index.indices.at(n)
							){
								return false;
							}
						}
					}
					else if(
						index.indices.at(n).isLabeledWildcard()
					){
						for(
							unsigned int c = 0;
							c < indices.size();
							c++
						){
							if(
								index.indices.at(c)
								== index.indices.at(n)
								&& indices.at(c)
								!= indices.at(n)
							){
								return false;
							}
						}
					}
					else{
						return false;
					}
				}
			}
		}
	}
	else{
		return false;
	}

	return true;
}

inline Subindex& Index::at(unsigned int n){
	return indices.at(n);
}

inline const Subindex& Index::at(unsigned int n) const{
	return indices.at(n);
}

inline unsigned int Index::getSize() const{
	return indices.size();
}

inline void Index::reserve(unsigned int size){
	indices.reserve(size);
}

inline void Index::pushBack(Subindex subindex){
	indices.push_back(subindex);
}

inline Subindex Index::popFront(){
	Subindex first = indices.at(0);
	indices.erase(indices.begin());

	return first;
}

inline Subindex Index::popBack(){
	Subindex last = indices.back();
	indices.pop_back();

	return last;
}

inline void Index::insert(unsigned int n, Subindex subindex){
	indices.insert(indices.begin() + n, subindex);
}

inline Subindex Index::erase(unsigned int n){
	Subindex subindex = indices[n];
	indices.erase(indices.begin() + n);

	return subindex;
}

inline std::vector<Index> Index::split() const{
	std::vector<Index> components;
	components.push_back(Index());
	for(unsigned int n = 0; n < indices.size(); n++){
		if(indices[n].isIndexSeparator())
			components.push_back(Index());
		else
			components.back().pushBack(indices[n]);
	}

	return components;
}

inline bool Index::isPatternIndex() const{
	for(unsigned int n = 0; n < indices.size(); n++)
		if(indices.at(n) < 0)
			return true;

	return false;
}

inline Subindex& Index::operator[](unsigned int subindex){
	return indices[subindex];
}

inline const Subindex& Index::operator[](unsigned int subindex) const{
	return indices[subindex];
}

inline unsigned int Index::getSizeInBytes() const{
	return sizeof(*this) + sizeof(int)*indices.capacity();
}

};	//End of namespace MyTBTK

#endif
