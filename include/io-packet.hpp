#pragma once

#import <cereal/archives/portable_binary.hpp>
#import <cereal/types/array.hpp>
#import <cereal/types/string.hpp>
#import <cereal/types/unordered_map.hpp>
#import <cereal/types/vector.hpp>
#import <string>
#import <unordered_map>
#import <vector>

#import "params.hpp"

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