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

/** @file OverlapAmplitudeSet.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/OverlapAmplitudeSet.h"
#include "MyTBTK/Streams.h"
#include "MyTBTK/MyTBTKMacros.h"

#include "MyTBTK/json.hpp"

using namespace std;
//using namespace nlohmann;

namespace MyTBTK{

OverlapAmplitudeSet::OverlapAmplitudeSet(){
    assumeOrthonormalBasis = true;
}

OverlapAmplitudeSet::OverlapAmplitudeSet(
    const string &serialization,
    Mode mode
){
    switch(mode){
    case Mode::JSON:
    {
        try{
            nlohmann::json j = nlohmann::json::parse(serialization);
            overlapAmplitudeTree
                = IndexedDataTree<OverlapAmplitude>(
                    j.at("overlapAmplitudeTree").dump(),
                    mode
                );
            assumeOrthonormalBasis
                = j.at("assumeOrthonormalBasis").get<bool>();
        }
        catch(nlohmann::json::exception &e){
            MyTBTKExit(
                "OverlapAmplitudeSet::OverlapAmplitudeSet()",
                "Unable to parse string as OverlapAmplitudeSet"
                << " '" << serialization << "'.",
                ""
            );
        }

        break;
    }
    default:
        MyTBTKExit(
            "OverlapAmplitudeSet::OverlapAmplitudeSet()",
            "Only Serializable::Mode::JSON is supported yet.",
            ""
        );
    }
}

OverlapAmplitudeSet::~OverlapAmplitudeSet(){
}

OverlapAmplitudeSet::Iterator OverlapAmplitudeSet::begin(){
    return OverlapAmplitudeSet::Iterator(overlapAmplitudeTree);
}

OverlapAmplitudeSet::ConstIterator OverlapAmplitudeSet::begin() const{
    return OverlapAmplitudeSet::ConstIterator(overlapAmplitudeTree);
}

OverlapAmplitudeSet::ConstIterator OverlapAmplitudeSet::cbegin() const{
    return OverlapAmplitudeSet::ConstIterator(overlapAmplitudeTree);
}

OverlapAmplitudeSet::Iterator OverlapAmplitudeSet::end(){
    return OverlapAmplitudeSet::Iterator(overlapAmplitudeTree, true);
}

OverlapAmplitudeSet::ConstIterator OverlapAmplitudeSet::end() const{
    return OverlapAmplitudeSet::ConstIterator(overlapAmplitudeTree, true);
}

OverlapAmplitudeSet::ConstIterator OverlapAmplitudeSet::cend() const{
    return OverlapAmplitudeSet::ConstIterator(overlapAmplitudeTree, true);
}

string OverlapAmplitudeSet::serialize(Mode mode) const{
    switch(mode){
    case Mode::JSON:
    {
        nlohmann::json j;
        j["id"] = "OverlapAmplitudeSet";
        j["overlapAmplitudeTree"] = nlohmann::json::parse(
            overlapAmplitudeTree.serialize(mode)
        );
        j["assumeOrthonormalBasis"]
            = nlohmann::json(assumeOrthonormalBasis);

        return j.dump();
    }
    default:
        MyTBTKExit(
            "OverlapAmplitudeSet::serialize()",
            "Only Serializable::Mode::JSON is supported yet.",
            ""
        );
    }
}

};    //End of namespace MyTBTK
