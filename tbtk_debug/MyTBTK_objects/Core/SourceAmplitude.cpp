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

/** @file SourceAmplitude.cpp
 *
 *  @author Kristofer Björnson
 */

#include "MyTBTK/SourceAmplitude.h"

#include "MyTBTK/json.hpp"

using namespace std;
//using namespace nlohmann;

namespace MyTBTK{

SourceAmplitude::SourceAmplitude(
    complex<double> amplitude,
    Index index
) :
    amplitude(amplitude),
    amplitudeCallback(nullptr),
    index(index)
{
};

SourceAmplitude::SourceAmplitude(
    const AmplitudeCallback &amplitudeCallback,
    Index index
) :
    amplitudeCallback(&amplitudeCallback),
    index(index)
{
};

SourceAmplitude::SourceAmplitude(
    const string &serialization,
    Serializable::Mode mode
){
    MyTBTKAssert(
        Serializable::validate(
            serialization,
            "SourceAmplitude",
            mode
        ),
        "SourceAmplitude::SourceAmplitude()",
        "Unable to parse string as SourceAmplitude '" << serialization
        << "'.",
        ""
    );

    switch(mode){
    case Serializable::Mode::JSON:
    {
        try{
            amplitudeCallback = nullptr;

            nlohmann::json j = nlohmann::json::parse(serialization);
            amplitude = Serializable::deserialize<complex<double>>(
                j["amplitude"].get<string>(),
                mode
            );
            index = Index(j["index"].dump(), mode);
        }
        catch(nlohmann::json::exception &e){
            MyTBTKExit(
                "SourceAmplitude::SourceAmplitude()",
                "Unable to parse string as SourceAmplitude '"
                << serialization << "'.",
                ""
            );
        }

        break;
    }
    default:
        MyTBTKExit(
            "SourceAmplitude::SourceAmplitude()",
            "Only Serializable::Mode::JSON is supported yet.",
            ""
        );
    }
}

string SourceAmplitude::serialize(Serializable::Mode mode) const{
    MyTBTKAssert(
        amplitudeCallback == nullptr,
        "SourceAmplitude::serialize()",
        "Unable to serialize SourceAmplitude that uses callback."
        << " value.",
        ""
    );

    switch(mode){
    case Serializable::Mode::JSON:
    {
        nlohmann::json j;
        j["id"] = "SourceAmplitude";
        j["amplitude"] = Serializable::serialize(amplitude, mode);
        j["index"] = nlohmann::json::parse(index.serialize(mode));

        return j.dump();
    }
    default:
        MyTBTKExit(
            "SourceAmplitude::serialize()",
            "Only Serializable::Mode::JSON is supported yet.",
            ""
        );
    }
}

};    //End of namespace MyTBTK
