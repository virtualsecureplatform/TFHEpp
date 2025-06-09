export module tfhepp:io_packet;
import std;

// Cereal imports - noted for special handling with modules
#import <cereal/archives/portable_binary.hpp>
#import <cereal/types/array.hpp>
#import <cereal/types/string.hpp>
#import <cereal/types/unordered_map.hpp>
#import <cereal/types/vector.hpp>

export import tfhepp:params;
export import tfhepp:tlwe;
export import tfhepp:trlwe;
export import tfhepp:trgsw;

export namespace TFHEpp {

export struct IOpacket {
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
        archive(tlwelvl0ports, tlwelvl1ports, trlwelvl2ports, trgswlvl2ports, numCycles); // Added tlwelvl1ports to serialization
    }
};
}  // namespace TFHEpp