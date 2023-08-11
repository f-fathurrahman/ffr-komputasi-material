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

/** @package MyTBTKcalc
 *  @file CArray.h
 *  @brief Container for a C style array.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_C_ARRAY
#define COM_DAFER45_MyTBTK_C_ARRAY

#include "MyTBTK/Serializable.h"
#include "MyTBTK/MyTBTKMacros.h"

#include "MyTBTK/json.hpp"

namespace MyTBTK{

/** @brief Container for a C style array.
 *
 *  The CArray contains a pointer to an array and its size. It provides the
 *  same efficiency as a raw array, but without requiring manual memory
 *  management.
 *
 *  # Example
 *  \snippet Utilities/CArray.cpp CArray
 *  ## Output
 *  \snippet output/Utilities/CArray.txt CArray */
template<typename DataType>
class CArray{
public:
	//MyTBTKFeature Utilities.CArray.construction.1 2019-10-30
	/** Constructor. */
	CArray();

	//MyTBTKFeature Utilities.CArray.construction.2 2019-10-30
	/** Constructor.
	 *
	 *  @param size The size of the array. */
	CArray(unsigned int size);

	//MyTBTKFeature Utilities.CArray.copyConstruction.1.C++ 2019-10-30
	/** Copy constructor.
	 *
	 *  @param carray The carray to copy. */
	CArray(const CArray &carray);

	//MyTBTKFeature Utilities.CArray.moveConstruction.1.C++ 2019-10-30
	/** Move constructor.
	 *
	 *  @param carray The carray to move. */
	CArray(CArray &&carray);

	/** Constructs a CArray from a serialization string.
	 *
	 *  @param serialization Serialization string from which to construct
	 *  the CArray.
	 *
	 *  @param mode Mode with which the string has been serialized. */
	CArray(const std::string &serialization, Serializable::Mode mode);

	/** Destructor. */
	~CArray();

	/** Assignment operator.
	 *
	 *  @param rhs The right hand side of the expression.
	 *
	 *  @return The left hand side after assignment has occured. */
	CArray& operator=(const CArray &carray);

	/** Move assignment operator.
	 *
	 *  @param rhs The right hand side of the expression.
	 *
	 *  @return The left hand side after assignment has occured. */
	CArray& operator=(CArray &&carray);

	//MyTBTKFeature Utilities.CArray.operatorArraySubscript.1 2019-10-30
	/** Array subscript operator.
	 *
	 *  @param n Position to get the value for.
	 *
	 *  @return The value at position n. */
	DataType& operator[](unsigned int n);

	//MyTBTKFeature Utilities.CArray.operatorArraySubscript.2 2019-10-30
	/** Array subscript operator.
	 *
	 *  @param n Position to get the value for.
	 *
	 *  @return The value at position n. */
	const DataType& operator[](unsigned int n) const;

	//MyTBTKFeature Utilities.CArray.getData.1.C++ 2019-10-31
	/** Get the data as a bare c-array.
	 *
	 *  @return Pointer to the data array. */
	DataType* getData();

	//MyTBTKFeature Utilities.CArray.getData.2.C++ 2019-10-31
	/** Get the data as a bare c-array.
	 *
	 *  @return Pointer to the data array. */
	const DataType* getData() const;

	/** Get size.
	 *
	 *  @return The size of the array. */
	unsigned int getSize() const;

	/** Serialize the CArray. Note that CArray is pseudo-Serializable in
	 *  that it implements Serializable interface, but does so
	 *  non-virtually.
	 *
	 *  @param mode Serialization mode to use.
	 *
	 *  @return Serialized string representation of the CArray. */
	std::string serialize(Serializable::Mode mode) const;
private:
	/** Size. */
	unsigned int size;

	/** Data. */
	DataType *data;
};

template<typename DataType>
CArray<DataType>::CArray(){
	data = nullptr;
}

template<typename DataType>
CArray<DataType>::CArray(unsigned int size){
	this->size = size;
	data = new DataType[size];
}

template<typename DataType>
CArray<DataType>::CArray(const CArray &carray){
	size = carray.size;
	if(carray.data == nullptr){
		data = nullptr;
	}
	else{
		data = new DataType[size];
		for(unsigned int n = 0; n < size; n++)
			data[n] = carray.data[n];
	}
}

template<typename DataType>
CArray<DataType>::CArray(CArray &&carray){
	size = carray.size;
	if(carray.data == nullptr){
		data = nullptr;
	}
	else{
		data = carray.data;
		carray.data = nullptr;
	}
}

template<typename DataType>
CArray<DataType>::CArray(
	const std::string &serialization,
	Serializable::Mode mode
){
	MyTBTKAssert(
		Serializable::validate(
			serialization,
			"CArray",
			mode
		),
		"CArray::CArray()",
		"Unable to parse string as CArray '" << serialization << "'.",
		""
	);

	switch(mode){
	case Serializable::Mode::JSON:
	{
		try{
			nlohmann::json j
				= nlohmann::json::parse(serialization);
			size = j.at("size").get<unsigned int>();
			nlohmann::json d = j.at("data");
			std::vector<DataType> tempData;
			for(
				nlohmann::json::iterator iterator = d.begin();
				iterator != d.end();
				++iterator
			){
				tempData.push_back(
					Serializable::deserialize<DataType>(
						*iterator,
						mode
					)
				);
			}
			MyTBTKAssert(
				size == tempData.size(),
				"CArray::CArray()",
				"Unable to deserialize CArray. The number of"
				<< "data elements does not agree with the size"
				<< " '" << serialization << "'.",
				""
			);
			data = new DataType[size];
			for(unsigned int n = 0; n < size; n++)
				data[n] = tempData[n];
		}
		catch(nlohmann::json::exception &e){
			MyTBTKExit(
				"CArray::CArray()",
				"Unable to parse string as CArray '"
				<< serialization << "'.",
				""
			);
		}

		break;
	}
	default:
		MyTBTKExit(
			"CArray::CArray()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

template<typename DataType>
CArray<DataType>::~CArray(){
	if(data != nullptr)
		delete [] data;
}

template<typename DataType>
CArray<DataType>& CArray<DataType>::operator=(const CArray &rhs){
	if(this != &rhs){
		size = rhs.size;
		if(data != nullptr)
			delete [] data;

		if(rhs.data == nullptr){
			data = nullptr;
		}
		else{
			data = new DataType[size];
			for(unsigned int n = 0; n < size; n++)
				data[n] = rhs.data[n];
		}
	}

	return *this;
}

template<typename DataType>
CArray<DataType>& CArray<DataType>::operator=(CArray &&rhs){
	if(this != &rhs){
		size = rhs.size;
		if(data != nullptr)
			delete [] data;

		if(rhs.data == nullptr){
			data = nullptr;
		}
		else{
			data = rhs.data;
			rhs.data = nullptr;
		}
	}

	return *this;
}

template<typename DataType>
DataType& CArray<DataType>::operator[](unsigned int n){
	return data[n];
}

template<typename DataType>
const DataType& CArray<DataType>::operator[](unsigned int n) const{
	return data[n];
}

template<typename DataType>
DataType* CArray<DataType>::getData(){
	return data;
}

template<typename DataType>
const DataType* CArray<DataType>::getData() const{
	return data;
}

template<typename DataType>
unsigned int CArray<DataType>::getSize() const{
	return size;
}

template<typename DataType>
std::string CArray<DataType>::serialize(Serializable::Mode mode) const{
	switch(mode){
	case Serializable::Mode::JSON:
	{
		nlohmann::json j;
		j["id"] = "CArray";
		j["size"] = size;
		j["data"] = nlohmann::json();
		for(unsigned int n = 0; n < size; n++){
			j["data"].push_back(
				Serializable::serialize(data[n], mode)
			);
		}

		return j.dump();
	}
	default:
		MyTBTKExit(
			"CArray::serialize()",
			"Only Serializable::Mode::JSON is supported yet.",
			""
		);
	}
}

}; //End of namesapce MyTBTK

#endif
