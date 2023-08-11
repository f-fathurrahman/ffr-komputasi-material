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
 *  @file Range.h
 *  @brief Helper class for flattening nested looping.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_MULTI_COUNTER
#define COM_DAFER45_MyTBTK_MULTI_COUNTER

#include "MyTBTK/MyTBTKMacros.h"

#include <vector>

namespace MyTBTK{

/** @brief Helper class for flattening nested loops.
 *
 *  The MultiCounter allows for multiple loop variables to be looped over using
 *  a single loop. It can be used to flatten deeply nested loops.
 *
 *  # Example
 *  \snippet Utilities/MultiCounter.cpp MultiCounter
 *  ## Output
 *  \snippet output/Utilities/MultiCounter.txt MultiCounter */
template<typename DataType>
class MultiCounter{
public:
	/** Constructor.
	 *
	 *  @param begin The start values for the counters.
	 *  @param end The end of the counters. This value is not included in
	 *  the range.
	 *
	 *  @param increment The increment size of the counter. */
	MultiCounter(
		const std::vector<DataType> &begin,
		const std::vector<DataType> &end,
		const std::vector<DataType> &increment
	);

	/** Increment operator.
	 *
	 *  @return A reference to the MultiCounter after the increment has
	 *  occured. */
	MultiCounter& operator++();

	/** Array subscript operator. Returns the current valu of the nth
	 *  counter.
	 *
	 *  @param n The counter to get the value for.
	 *
	 *  @return The value of the nth counter. */
	const DataType operator[](unsigned int n) const;

	/** */
	operator std::vector<DataType>() const;

	/** Reset the counter. */
	void reset();

	/** Check if the counter has reached the end.
	 *
	 *  @return True if the counter has reached the end, otherwise false.
	 */
	bool done() const;

	/** Get the number of counters.
	 *
	 *  @return The number of counters. */
	unsigned int getSize() const;
private:
	/** Values at which the iteration begins. */
	std::vector<DataType> begin;

	/** Values at which the iteration ends. */
	std::vector<DataType> end;

	/** Increments for the iteration. */
	std::vector<DataType> increment;

	/** Counters for the iteration. */
	std::vector<DataType> counter;
};

template<typename DataType>
inline MultiCounter<DataType>::MultiCounter(
	const std::vector<DataType> &begin,
	const std::vector<DataType> &end,
	const std::vector<DataType> &increment
){
	MyTBTKAssert(
		begin.size() == end.size(),
		"MultiCounter::MultiCounter()",
		"'begin' and 'end' must have the same number of elements.",
		""
	);
	MyTBTKAssert(
		begin.size() == increment.size(),
		"MultiCounter::MultiCounter()",
		"'begin' and 'increment' must have the same number of"
		<< " elements.",
		""
	);

	this->begin = begin;
	this->end = end;
	this->increment = increment;

	for(unsigned int n = 0; n < begin.size(); n++){
		MyTBTKAssert(
			this->begin[n] <= this->end[n],
			"MultiCounter::MultiCounter()",
			"Only forward iteration is supported, but entry '" << n
			<< "' in 'end' is smaller than entry '" << n << "' in"
			<< " 'begin'.",
			""
		);

		MyTBTKAssert(
			this->increment[n] > 0,
			"MultiCounter::MultiCounter()",
			"Only positive increments are supported, but entry '"
			<< n << "' in 'increment' is '" << this->increment[n]
			<< "'.",
			""
		);
	}

	reset();
}

template<typename DataType>
inline MultiCounter<DataType>& MultiCounter<DataType>::operator++(){
	for(int n = counter.size()-1; n > 0; n--){
		counter[n] += increment[n];
		if(counter[n] < end[n])
			return *this;

		counter[n] = begin[n];
	}

	counter[0] += increment[0];

	return *this;
}

template<typename DataType>
inline const DataType MultiCounter<DataType>::operator[](unsigned int n) const{
	return counter[n];
}

template<typename DataType>
inline MultiCounter<DataType>::operator std::vector<DataType>() const{
	return counter;
}

template<typename DataType>
inline void MultiCounter<DataType>::reset(){
	counter = begin;
}

template<typename DataType>
inline bool MultiCounter<DataType>::done() const{
	return counter[0] >= end[0];
}

template<typename DataType>
inline unsigned int MultiCounter<DataType>::getSize() const{
	return begin.size();
}

}; //End of namesapce MyTBTK

#endif
