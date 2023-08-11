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

/** @package MyTBTKcalc
 *  @file IndexDescriptor.h
 *  @brief Describes the index structure of data stored for several indices.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_INDEX_DESCRIPTOR
#define COM_DAFER45_MyTBTK_INDEX_DESCRIPTOR

#include "MyTBTK/IndexedDataTree.h"
#include "MyTBTK/IndexException.h"
#include "MyTBTK/IndexTree.h"
#include "MyTBTK/Serializable.h"
#include "MyTBTK/MyTBTKMacros.h"

namespace MyTBTK{

/** @brief Describes the index structure of data stored for several indices.
 *
 *  The IndexDescriptor is a helper class for the AbstractProperty to help
 *  handle the storage of data using different type of indexing. See
 *  AbstractProperty for a description of the various formats. */
class IndexDescriptor : public Serializable{
public:
	/** Enum class determining the storage format. */
	enum class Format {None, Ranges, Custom, Dynamic};

	/** Constructs an IndexDescriptor on the Format::None. */
	IndexDescriptor();

	/** Constructs an IndexDescriptor on the Format::Ranges.
	 *
	 *  @param ranges The upper limits (exlcusive) of the grid described by
	 *  the IndexDescriptor. */
	IndexDescriptor(const std::vector<int> &ranges);

	/** Constructs an IndexDescriptor on the Format::Custom.
	 *
	 *  @param indexTree An IndexTree containing all the @link Index
	 *  Indices @endlink that should be described by the IndexDescriptor.
	 */
	IndexDescriptor(const IndexTree &indexTree);

	/** Copy constructor.
	 *
	 *  @param indexDescriptor IndexDescriptor to copy. */
	IndexDescriptor(const IndexDescriptor &indexDescriptor);

	/** Move constructor.
	 *
	 *  @param indexDescriptor IndexDescriptor to move. */
	IndexDescriptor(IndexDescriptor &&indexDescriptor);

	/** Constructor. Construct the IndexDescriptor from a serialization
	 *  string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the IndexDescriptor.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	IndexDescriptor(const std::string &serialization, Mode mode);

	/** Destructor. */
	~IndexDescriptor();

	/** Assignment operator.
	 *
	 *  @param rhs IndexDescriptor to assign to the left hand side.
	 *
	 *  @return Reference to the assigned IndexDescriptor. */
	IndexDescriptor& operator=(const IndexDescriptor &rhs);

	/** Move assignment operator.
	 *
	 *  @param rhs IndexDescriptor to assign to the left hand side.
	 *
	 *  @return Reference to the assigned IndexDescriptor. */
	IndexDescriptor& operator=(IndexDescriptor &&rhs);

	/** Comparison operator.
	 *
	 *  @param rhs Left hand side of the equality sign.
	 *  @param rhs Right hand side of the equality sign.
	 *
	 *  @return True if the @link IndexDescriptor IndexDescriptors @endlink
	 *  describe identical index structures, otherwise false. */
	friend bool operator==(
		const IndexDescriptor &lhs,
		const IndexDescriptor &rhs
	);

	/** Inequality operator.
	 *
	 *  @param rhs Left hand side of the inequality sign.
	 *  @param rhs Right hand side of the inequality sign.
	 *
	 *  @return False if the @link IndexDescriptor IndexDescriptors @endlink
	 *  describe identical index structures, otherwise true. */
	friend bool operator!=(
		const IndexDescriptor &lhs,
		const IndexDescriptor &rhs
	);

	/** Get format.
	 *
	 *  @return The Format that the IndexDescriptor describes. */
	Format getFormat() const;

	/** Get ranges. [Only works for the Ranges format.]
	 *
	 *  @return The upper limits (exclusive) of the grid described by the
	 *  IndexDescriptor. */
	std::vector<int> getRanges() const;

	/** Get IndexTree. [Only works for the Custom format.]
	 *
	 *  @return An IndexTree containing all the @link Index Indices
	 *  @endlink that are described by the IndexDescriptor. */
	const IndexTree& getIndexTree() const;

	/** Add Index. [Only works for the Dynamic format.]
	 *
	 *  @param index Index to add.  */
	void add(const Index &index);

	/** Get linear index. [Only works for the Custom and Dynamic formats.]
	 *
	 *  @param index The Index to get the linear Index for.
	 *  @param returnNegativeForMissingIndex If set to true, requesting an
	 *  Index that is not included results in -1 being returned. If set to
	 *  false, an IndexException will be thrown in the same situation.
	 *
	 *  @return The linear index that corresponds to the given Index. */
	int getLinearIndex(
		const Index &index,
		bool returnNegativeForMissingIndex = false
	) const;

	/** Get the number of data elements described by the IndexDescriptor.
	 *  For Format::None this is equal to 1, while for the other formats
	 *  this is equal to the number of @link Index Indices @endlink that
	 *  are described by the IndexDescriptor.
	 *
	 *  @return The number of data elements described by the
	 *  IndexDescriptor. */
	unsigned int getSize() const;

	/** Check whether a given Index is contained in the IndexDescriptor.
	 *  [Only works for the Custom and Dynamic formats.]
	 *
	 *  @return True if the index descriptor contains the given index. */
	bool contains(const Index &index) const;

	/** Implements Serializable::serialize(). */
	virtual std::string serialize(Mode mode) const;
private:
	/** Index descriptor format. */
	Format format;

	class NoneFormat{
	public:
	};

	class RangeFormat{
	public:
		/** Number of dimensions. */
		unsigned int dimensions;

		/** Ranges. */
		int *ranges;
	};

	class CustomFormat{
	public:
		/** IndexTree. */
		IndexTree *indexTree;
	};

	class DynamicFormat{
	public:
		/** IndexedDataTree that maps each added Index to a linear
		 *  index. */
		IndexedDataTree<unsigned int> *indexedDataTree;

		/** Current number of Indices added to the IndexedDataTree. */
		unsigned int size;
	};

	/** Union of descriptor formats. */
	union Descriptor{
		NoneFormat noneFormat;
		RangeFormat rangeFormat;
		CustomFormat customFormat;
		DynamicFormat dynamicFormat;
	};

	/** Actuall descriptor. */
	Descriptor descriptor;
};

inline bool operator!=(const IndexDescriptor &lhs, const IndexDescriptor &rhs){
	return !(rhs == lhs);
}

inline IndexDescriptor::Format IndexDescriptor::getFormat() const{
	return format;
}

inline std::vector<int> IndexDescriptor::getRanges() const{
	MyTBTKAssert(
		format == Format::Ranges,
		"IndexDescriptor::setDimensions()",
		"The IndexDescriptor is not of the format Format::Ranges.",
		""
	);

	std::vector<int> ranges;
	for(unsigned int n = 0; n < descriptor.rangeFormat.dimensions; n++)
		ranges.push_back(descriptor.rangeFormat.ranges[n]);

	return ranges;
}

inline const IndexTree& IndexDescriptor::getIndexTree() const{
	MyTBTKAssert(
		format == Format::Custom,
		"IndexDescriptor::getIndexTree()",
		"The IndexDescriptor is not of the format Format::Custom.",
		""
	);

	return *descriptor.customFormat.indexTree;
}

inline void IndexDescriptor::add(const Index &index){
	MyTBTKAssert(
		format == Format::Dynamic,
		"IndexDescriptor::add()",
		"The IndexDescriptor is not of the format Format::Dynamic.",
		""
	);

	unsigned int dummy;
	MyTBTKAssert(
		!descriptor.dynamicFormat.indexedDataTree->get(dummy, index),
		"IndexDescriptor::add()",
		"The IndexDescriptor already contains the Index '"
		<< index.toString() << "'.",
		""
	);

	descriptor.dynamicFormat.indexedDataTree->add(
		descriptor.dynamicFormat.size,
		index
	);
	descriptor.dynamicFormat.size++;
}

inline int IndexDescriptor::getLinearIndex(
	const Index &index,
	bool returnNegativeForMissingIndex
) const{
	switch(format){
	case Format::Custom:
		return descriptor.customFormat.indexTree->getLinearIndex(
			index,
			IndexTree::SearchMode::MatchWildcards,
			returnNegativeForMissingIndex
		);
	case Format::Dynamic:
	{
		unsigned int linearIndex;
		if(
			descriptor.dynamicFormat.indexedDataTree->get(
				linearIndex,
				index
			)
		){
			return linearIndex;
		}
		else if(returnNegativeForMissingIndex){
			return -1;
		}
		else{
			throw IndexException(
				"IndexDescriptor::getLinearIndex()",
				MyTBTKWhere,
				"Index not included in the IndexDescriptor '"
				+ index.toString() + "'.",
				""
			);
		}
	}
	default:
		MyTBTKExit(
			"IndexDescriptor::getOffset()",
			"The IndexDescriptor is not of the format"
			" Format::Custom or Format::Dynamic.",
			""
		);
	}
}

inline bool IndexDescriptor::contains(const Index &index) const{
	switch(format){
	case Format::Custom:
		if(
			descriptor.customFormat.indexTree->getLinearIndex(
				index,
				IndexTree::SearchMode::StrictMatch,
				true
			) == -1
		){
			return false;
		}
		else{
			return true;
		}
	case Format::Dynamic:
	{
		unsigned int dummy;
		if(
			descriptor.dynamicFormat.indexedDataTree->get(
				dummy,
				index
			)
		){
			return true;
		}
		else{
			return false;
		}
	}
	default:
		MyTBTKExit(
			"IndexDescriptor::contains()",
			"The IndexDescriptor is not of the format"
			" Format::Custom or Format::Dynamic.",
			""
		);
	}
}

};	//End namespace MyTBTK

#endif
