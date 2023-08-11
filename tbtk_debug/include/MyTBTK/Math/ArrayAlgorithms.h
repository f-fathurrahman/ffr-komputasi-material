/* Copyright 2020 Kristofer Björnson
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
 *  @file ArrayAlgorithms.h
 *  @brief Implements Array algorithms.
 *
 *  @author Kristofer Björnson
 */

#ifndef COM_DAFER45_MyTBTK_MATH_ARRAY_ALGORITHMS
#define COM_DAFER45_MyTBTK_MATH_ARRAY_ALGORITHMS

#include "MyTBTK/Array.h"

#include <cmath>

namespace MyTBTK{
namespace Math{

template<typename DataType>
class ArrayAlgorithms{
public:
	/** Elementwise sine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise sine of the input Array. */
	static Array<DataType> sin(const Array<DataType> &array);

	/** Elementwise cosine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise cosine of the input Array. */
	static Array<DataType> cos(const Array<DataType> &array);

	/** Elementwise tangens.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise tangens of the input Array. */
	static Array<DataType> tan(const Array<DataType> &array);

	/** Elementwise arcsine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise arcsine of the input Array. */
	static Array<DataType> asin(const Array<DataType> &array);

	/** Elementwise arccosine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise arccosine of the input Array. */
	static Array<DataType> acos(const Array<DataType> &array);

	/** Elementwise arctangens.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise arctangens of the input Array. */
	static Array<DataType> atan(const Array<DataType> &array);

	/** Elementwise hyperbolic sine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise hyperbolic sine of the input Array. */
	static Array<DataType> sinh(const Array<DataType> &array);

	/** Elementwise hyperbolic cosine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise hyperbolic cosine of the input Array. */
	static Array<DataType> cosh(const Array<DataType> &array);

	/** Elementwise hyperbolic tangens.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise hyperbolic tangens of the input Array. */
	static Array<DataType> tanh(const Array<DataType> &array);

	/** Elementwise hyperbolic arcsine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise hyperbolic arcsine of the input Array. */
	static Array<DataType> asinh(const Array<DataType> &array);

	/** Elementwise hyperbolic arccosine.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise hyperbolic arccosine of the input Array. */
	static Array<DataType> acosh(const Array<DataType> &array);

	/** Elementwise hyperbolic arctangens.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise hyperbolic arctangens of the input Array. */
	static Array<DataType> atanh(const Array<DataType> &array);

	/** Elementwise natural logarithm.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise natural logarithm of the input Array. */
	static Array<DataType> log(const Array<DataType> &array);

	/** Elementwise base-2 logarithm.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise base-2 logarithm of the input Array. */
	static Array<DataType> log2(const Array<DataType> &array);

	/** Elementwise base-10 logarithm.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise base-10 logarithm of the input Array. */
	static Array<DataType> log10(const Array<DataType> &array);

	/** Elementwise exponent.
	 *
	 *  @param array Input Array.
	 *  @param exponent The exponent to rise the elements to.
	 *
	 *  @return The elementwise exponent of the input Array. */
	static Array<DataType> pow(
		const Array<DataType> &array,
		double exponent
	);

	/** Elementwise exponential.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise exponential of the input Array. */
	static Array<DataType> exp(const Array<DataType> &array);

	/** Elementwise absolute value.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise absolute value of the input Array. */
	template<typename T = DataType>
	static typename std::enable_if<
		!std::is_same<T, std::complex<double>>::value,
		Array<DataType>
	>::type abs(const Array<DataType> &array);

	/** Elementwise absolute value.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise absolute value of the input Array. */
	template<typename T = DataType>
	static typename std::enable_if<
		std::is_same<T, std::complex<double>>::value,
		Array<double>
	>::type abs(const Array<DataType> &array);

	/** Elementwise argument of complex number.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise argument of complex number of the input
	 *  Array. */
	template<typename T = DataType>
	static typename std::enable_if<
		std::is_same<T, std::complex<double>>::value,
		Array<double>
	>::type arg(const Array<std::complex<double>> &array);

	/** Elementwise real component of complex number.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise real component of complex value of the input
	 *  Array. */
	template<typename T = DataType>
	static typename std::enable_if<
		std::is_same<T, std::complex<double>>::value,
		Array<double>
	>::type real(const Array<std::complex<double>> &array);

	/** Elementwise imaginary component of complex number.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise imaginary component of complex number of the
	 *  input Array. */
	template<typename T = DataType>
	static typename std::enable_if<
		std::is_same<T, std::complex<double>>::value,
		Array<double>
	>::type imag(const Array<std::complex<double>> &array);

	/** Elementwise complex conugate.
	 *
	 *  @param array Input Array.
	 *
	 *  @return The elementwise complex conjugate of the input Array. */
	template<typename T = DataType>
	static typename std::enable_if<
		std::is_same<T, std::complex<double>>::value,
		Array<DataType>
	>::type conj(const Array<DataType> &array);

	/** Elementwise square root.
	 *
	 * @param Input Array.
	 *
	 *  @return The elementwise square root of the input Array. */
	static Array<DataType> sqrt(const Array<DataType> &array);

	/** Maximum value.
	 *
	 * @param Input Array.
	 *
	 *  @return The maximum element of the input Array. */
	static DataType max(const Array<DataType> &array);

	/** Minimum value.
	 *
	 * @param Input Array.
	 *
	 *  @return The minimum element of the input Array. */
	static DataType min(const Array<DataType> &array);
private:
};

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::sin(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::sin(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::cos(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::cos(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::tan(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::tan(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::asin(
	const Array<DataType> &array
){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::asin(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::acos(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::acos(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::atan(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::atan(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::sinh(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::sinh(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::cosh(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::cosh(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::tanh(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::tanh(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::asinh(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::asinh(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::acosh(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::acosh(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::atanh(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::atanh(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::log(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::log(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::log2(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::log2(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::log10(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::log10(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::pow(
	const Array<DataType> &array,
	double exponent
){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::pow(arrayData[n], exponent);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::exp(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::exp(arrayData[n]);

	return result;
}

template<typename DataType>
template<typename T>
typename std::enable_if<
	!std::is_same<T, std::complex<double>>::value,
	Array<DataType>
>::type ArrayAlgorithms<DataType>::abs(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::abs(arrayData[n]);

	return result;
}

template<typename DataType>
template<typename T>
typename std::enable_if<
	std::is_same<T, std::complex<double>>::value,
	Array<double>
>::type ArrayAlgorithms<DataType>::abs(const Array<DataType> &array){
	Array<double> result = Array<double>::create(array.getRanges());
	CArray<double> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::abs(arrayData[n]);

	return result;
}

template<typename DataType>
template<typename T>
typename std::enable_if<
	std::is_same<T, std::complex<double>>::value,
	Array<double>
>::type ArrayAlgorithms<DataType>::arg(const Array<std::complex<double>> &array){
	Array<double> result = Array<double>::create(array.getRanges());
	CArray<double> &resultData = result.getData();
	const CArray<std::complex<double>> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::arg(arrayData[n]);

	return result;
}

template<typename DataType>
template<typename T>
typename std::enable_if<
	std::is_same<T, std::complex<double>>::value,
	Array<double>
>::type ArrayAlgorithms<DataType>::real(const Array<std::complex<double>> &array){
	Array<double> result = Array<double>::create(array.getRanges());
	CArray<double> &resultData = result.getData();
	const CArray<std::complex<double>> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::real(arrayData[n]);

	return result;
}

template<typename DataType>
template<typename T>
typename std::enable_if<
	std::is_same<T, std::complex<double>>::value,
	Array<double>
>::type ArrayAlgorithms<DataType>::imag(const Array<std::complex<double>> &array){
	Array<double> result = Array<double>::create(array.getRanges());
	CArray<double> &resultData = result.getData();
	const CArray<std::complex<double>> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::imag(arrayData[n]);

	return result;
}

template<typename DataType>
template<typename T>
typename std::enable_if<
	std::is_same<T, std::complex<double>>::value,
	Array<DataType>
>::type ArrayAlgorithms<DataType>::conj(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::conj(arrayData[n]);

	return result;
}

template<typename DataType>
Array<DataType> ArrayAlgorithms<DataType>::sqrt(const Array<DataType> &array){
	Array<DataType> result = Array<DataType>::create(array.getRanges());
	CArray<DataType> &resultData = result.getData();
	const CArray<DataType> &arrayData = array.getData();

	for(unsigned int n = 0; n < array.getSize(); n++)
		resultData[n] = std::sqrt(arrayData[n]);

	return result;
}

template<typename DataType>
DataType ArrayAlgorithms<DataType>::max(const Array<DataType> &array){
	const CArray<DataType> &arrayData = array.getData();
	DataType maximum = arrayData[0];
	for(unsigned int n = 1; n < array.getSize(); n++)
		if(maximum < arrayData[n])
			maximum = arrayData[n];

	return maximum;
}

template<typename DataType>
DataType ArrayAlgorithms<DataType>::min(const Array<DataType> &array){
	const CArray<DataType> &arrayData = array.getData();
	DataType minimum = arrayData[0];
	for(unsigned int n = 1; n < array.getSize(); n++)
		if(minimum > arrayData[n])
			minimum = arrayData[n];

	return minimum;
}

/** Elementwise sine.
 *
 * @param Input Array.
 *
 *  @return The elementwise sine of the input Array. */
template<typename DataType>
Array<DataType> sin(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::sin(array);
}

/** Elementwise cosine.
 *
 * @param Input Array.
 *
 *  @return The elementwise cosine of the input Array. */
template<typename DataType>
Array<DataType> cos(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::cos(array);
}

/** Elementwise tangens.
 *
 * @param Input Array.
 *
 *  @return The elementwise tangens of the input Array. */
template<typename DataType>
Array<DataType> tan(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::tan(array);
}

/** Elementwise arcsine.
 *
 * @param Input Array.
 *
 *  @return The elementwise arcsine of the input Array. */
template<typename DataType>
Array<DataType> asin(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::asin(array);
}

/** Elementwise arccosine.
 *
 * @param Input Array.
 *
 *  @return The elementwise arccosine of the input Array. */
template<typename DataType>
Array<DataType> acos(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::acos(array);
}

/** Elementwise arctangens.
 *
 * @param Input Array.
 *
 *  @return The elementwise arctangens of the input Array. */
template<typename DataType>
Array<DataType> atan(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::atan(array);
}

/** Elementwise hyperbolic sine.
 *
 * @param Input Array.
 *
 *  @return The elementwise hyperbolic sine of the input Array. */
template<typename DataType>
Array<DataType> sinh(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::sinh(array);
}

/** Elementwise hyperbolic cosine.
 *
 * @param Input Array.
 *
 *  @return The elementwise hyperbolic cosine of the input Array. */
template<typename DataType>
Array<DataType> cosh(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::cosh(array);
}

/** Elementwise hyperbolic tangens.
 *
 * @param Input Array.
 *
 *  @return The elementwise hyperbolic tangens of the input Array. */
template<typename DataType>
Array<DataType> tanh(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::tanh(array);
}

/** Elementwise hyperbolic arcsine.
 *
 * @param Input Array.
 *
 *  @return The elementwise hyperbolic arcsine of the input Array. */
template<typename DataType>
Array<DataType> asinh(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::asinh(array);
}

/** Elementwise hyperbolic arccosine.
 *
 * @param Input Array.
 *
 *  @return The elementwise hyperbolic arccosine of the input Array. */
template<typename DataType>
Array<DataType> acosh(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::acosh(array);
}

/** Elementwise hyperbolic arctangens.
 *
 * @param Input Array.
 *
 *  @return The elementwise hyperbolic arctangens of the input Array. */
template<typename DataType>
Array<DataType> atanh(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::atanh(array);
}

/** Elementwise natural logarithm.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise natural logarithm of the input Array. */
template<typename DataType>
Array<DataType> log(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::log(array);
}

/** Elementwise base-2 logarithm.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise base-2 logarithm of the input Array. */
template<typename DataType>
Array<DataType> log2(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::log2(array);
}

/** Elementwise base-10 logarithm.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise base-10 logarithm of the input Array. */
template<typename DataType>
Array<DataType> log10(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::log10(array);
}

/** Elementwise exponent.
 *
 *  @param array Input Array.
 *  @param exponent The exponent to rise the elements to.
 *
 *  @return The elementwise exponent of the input Array. */
template<typename DataType>
Array<DataType> pow(const Array<DataType> &array, double exponent){
	return ArrayAlgorithms<DataType>::pow(array, exponent);
}

/** Elementwise exponential.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise exponential of the input Array. */
template<typename DataType>
Array<DataType> exp(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::exp(array);
}

/** Elementwise absolute value.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise absolute value of the input Array. */
template<typename DataType>
typename std::enable_if<
	!std::is_same<DataType, std::complex<double>>::value,
	Array<DataType>
>::type abs(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::abs(array);
}

/** Elementwise absolute value.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise absolute value of the input Array. */
template<typename DataType>
typename std::enable_if<
	std::is_same<DataType, std::complex<double>>::value,
	Array<double>
>::type abs(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::abs(array);
}

/** Elementwise argument.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise argument of the input Array. */
inline Array<double> arg(const Array<std::complex<double>> &array){
	return ArrayAlgorithms<std::complex<double>>::arg(array);
}

/** Elementwise real component of complex number.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise real componenet of complex number of the input
 *  Array. */
inline Array<double> real(const Array<std::complex<double>> &array){
	return ArrayAlgorithms<std::complex<double>>::real(array);
}

/** Elementwise imaginary component of complex number.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise imaginary component of complex number of the input
 *  Array. */
inline Array<double> imag(const Array<std::complex<double>> &array){
	return ArrayAlgorithms<std::complex<double>>::imag(array);
}

/** Elementwise complex conjugate.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise complex conjugate of the input Array. */
inline Array<std::complex<double>> conj(const Array<std::complex<double>> &array){
	return ArrayAlgorithms<std::complex<double>>::conj(array);
}

/** Elementwise square root.
 *
 *  @param array Input Array.
 *
 *  @return The elementwise square root of the input Array. */
template<typename DataType>
Array<DataType> sqrt(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::sqrt(array);
}

/** Maximum value.
 *
 *  @param array Input Array.
 *
 *  @return The maximum element of the input Array. */
template<typename DataType>
DataType max(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::max(array);
}

/** Minimum value.
 *
 *  @param array Input Array.
 *
 *  @return The minimum element of the input Array. */
template<typename DataType>
DataType min(const Array<DataType> &array){
	return ArrayAlgorithms<DataType>::min(array);
}

}; //End of namespace Math
}; //End of namespace MyTBTK

#endif
