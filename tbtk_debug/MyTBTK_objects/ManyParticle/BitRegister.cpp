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

/** @file BitRegister.cpp
 *
 *  @author Kristofer Björnson
 */


#include "MyTBTK/BitRegister.h"
#include "MyTBTK/Streams.h"

namespace MyTBTK{

BitRegister::BitRegister(unsigned int numBits){
}

BitRegister::BitRegister(const BitRegister &bitRegister){
	values = bitRegister.values;
}

/*BitRegister::~BitRegister(){
}*/

};	//End of namespace MyTBTK