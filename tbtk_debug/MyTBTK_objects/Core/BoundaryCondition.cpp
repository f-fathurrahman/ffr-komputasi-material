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

/** @file BoundaryCondition.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/BoundaryCondition.h"

#include "MyTBTK/json.hpp"

using namespace std;
//using namespace nlohmann;

namespace MyTBTK{

BoundaryCondition::BoundaryCondition(
    const string &serialization,
    Serializable::Mode mode
){
    MyTBTKAssert(
        Serializable::validate(
            serialization,
            "BoundaryCondition",
            mode
        ),
        "BoundaryCondition::BoundaryCondition()",
        "Unable to parse string as BoundaryCondition '"
        << serialization << "'.",
        ""
    );

    switch(mode){
    case Serializable::Mode::JSON:
    {
        try{
            nlohmann::json j = nlohmann::json::parse(serialization);
            hoppingAmplitudeList = HoppingAmplitudeList(
                j["hoppingAmplitudeList"].dump(),
                mode
            );
            sourceAmplitude = SourceAmplitude(
                j["sourceAmplitude"].dump(),
                mode
            );
            eliminationIndex = Index(
                j["eliminationIndex"].dump(),
                mode
            );
        }
        catch(nlohmann::json::exception &e){
            MyTBTKExit(
                "BoundaryCondition::BoundaryCondition()",
                "Unable to parse string as BoundaryCondition '"
                << serialization << "'.",
                ""
            );
        }

        break;
    }
    default:
        MyTBTKExit(
            "BoundaryCondition::BoundaryCondition()",
            "Only Serializable::Mode::JSON is supported yet.",
            ""
        );
    }
}

string BoundaryCondition::serialize(Serializable::Mode mode) const{
    switch(mode){
    case Serializable::Mode::JSON:
    {
        nlohmann::json j;
        j["id"] = "BoundaryCondition";
        j["hoppingAmplitudeList"] = nlohmann::json::parse(hoppingAmplitudeList.serialize(mode));
        j["sourceAmplitude"] = nlohmann::json::parse(sourceAmplitude.serialize(mode));
        j["eliminationIndex"] = nlohmann::json::parse(eliminationIndex.serialize(mode));

        return j.dump();
    }
    default:
        MyTBTKExit(
            "BoundaryCondition::serialize()",
            "Only Serializable::Mode::JSON is supported yet.",
            ""
        );
    }
}

};    //End of namespace MyTBTK
