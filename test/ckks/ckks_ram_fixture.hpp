#pragma once
#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/tlwe_coefficient_packing.hpp>
#include <tfhe/refreshable_ram.hpp>

namespace ckks_ram {
using P=narrow_benchmark::BenchmarkParam<14,7>;
using KS=TFHEpp::independentlvl10param;
using CB=TFHEpp::cblvl02param;
using PKS=TFHEpp::cblvl21param;
using BR=TFHEpp::independentlvl01param;
using L=typename PKS::targetP;
using PackKey=TFHEpp::TLWECoefficientPackingKey48<P,L>;
using OutKey=TFHEpp::CKKSToTFHEKeySwitchKey<P,L,2,16>;
constexpr std::array<std::uint64_t,7> metadata{0x434b52414d303031ULL,P::n,448,48,L::n,2,16};
template<class... T>void save(const std::filesystem::path &file,const T&... values)
{
    std::ofstream f(file,std::ios::binary);if(!f)throw std::runtime_error("cannot create "+file.string());
    cereal::PortableBinaryOutputArchive ar(f);ar(metadata,values...);
    if(!f)throw std::runtime_error("fixture write failed");
}
template<class... T>void load(const std::filesystem::path &file,T&... values)
{
    std::ifstream f(file,std::ios::binary);if(!f)throw std::runtime_error("missing "+file.string());
    cereal::PortableBinaryInputArchive ar(f);std::array<std::uint64_t,7> meta;ar(meta);
    if(meta!=metadata)throw std::runtime_error("RAM fixture profile mismatch");
    ar(values...);if(f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("RAM fixture trailing data");
}
struct PlainRequest {std::uint32_t address;bool enable;std::uint32_t data;};
// Deliberately small first RAM: every address, disabled writes, repeated writes,
// and alternating bits. Word width and address depth are compile-time choices.
template<std::size_t A,std::size_t W>struct Fixture {
    using RAM=TFHEpp::RefreshableRAM<KS,CB,PKS,BR,A,W>;
    typename RAM::Memory initial;
    std::vector<typename RAM::Address> addresses;
    std::vector<TFHEpp::TLWE<L>> enables;
    std::vector<typename RAM::Word> data;
    template<class Archive>void serialize(Archive &ar){ar(initial,addresses,enables,data);}
};
inline std::vector<PlainRequest> trace()
{
    return {{0,true,1},{1,true,2},{2,true,3},{3,true,0},
            {0,false,2},{1,false,3},{2,false,0},{3,false,1},
            {0,true,2},{0,true,1},{1,true,1},{2,true,2},
            {0,false,0},{1,false,0},{2,false,0},{3,false,0}};
}
} // namespace ckks_ram
