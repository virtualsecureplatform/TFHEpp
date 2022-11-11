#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "params.hpp"

namespace TFHEpp {

struct IOpacket {
    std::unordered_map<std::string,
                       std::vector<TFHEpp::TLWE<TFHEpp::lvl0param>>>
        tlwelvl0ports;
    std::unordered_map<std::string,
                       std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>>>
        tlwelvl1ports;
    std::unordered_map<std::string,
                       std::vector<TFHEpp::TRLWE<TFHEpp::lvl2param>>>
        trlwelvl2ports;
    std::unordered_map<std::string,
                       std::vector<TFHEpp::TRGSW<TFHEpp::lvl2param>>>
        trgswlvl2ports;

    uint numCycles = 1;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(tlwelvl0ports, trlwelvl2ports, trgswlvl2ports, numCycles);
    }
};
}  // namespace TFHEpp