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
inline constexpr std::array<std::uint64_t,7> legacy_metadata{0x434b52414d303031ULL,P::n,448,48,L::n,2,16};
#ifdef USE_LVL1_SIGMA256
constexpr std::array<std::uint64_t,7> metadata{0x434b52414d303032ULL,P::n,448,48,L::n,2,16};
#else
constexpr auto metadata=legacy_metadata;
#endif
// SecretKey's default constructor generates secrets before deserialization.
// Read the identical archive layout without invoking that constructor.
struct SavedOwner {
    decltype(TFHEpp::lweKey::keys) keys{};
    TFHEpp::lweParams params;
    template<class Archive>void serialize(Archive &ar){ar(keys,params);}
    const TFHEpp::Key<L> &secret()const
    {return std::get<TFHEpp::KeyPair<TFHEpp::Key<L>>>(keys).independent;}
};
template<class... T>void save(const std::filesystem::path &file,const T&... values)
{
    std::ofstream f(file,std::ios::binary);if(!f)throw std::runtime_error("cannot create "+file.string());
    cereal::PortableBinaryOutputArchive ar(f);ar(metadata);
#ifdef USE_LVL1_SIGMA256
    ar(TFHEpp::TFHEArchiveNoiseProfile());
#endif
    ar(values...);
    if(!f)throw std::runtime_error("fixture write failed");
}
template<class... T>void load(const std::filesystem::path &file,T&... values)
{
    std::ifstream f(file,std::ios::binary);if(!f)throw std::runtime_error("missing "+file.string());
    cereal::PortableBinaryInputArchive ar(f);std::array<std::uint64_t,7> meta;ar(meta);
    if(meta!=metadata)throw std::runtime_error("RAM fixture profile mismatch");
#ifdef USE_LVL1_SIGMA256
    decltype(TFHEpp::TFHEArchiveNoiseProfile()) profile{};ar(profile);
    if(profile!=TFHEpp::TFHEArchiveNoiseProfile())throw std::runtime_error("RAM fixture noise mismatch");
#endif
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
template<std::size_t A=2,std::size_t W=2>
inline std::vector<PlainRequest> trace()
{
    if constexpr((A>=8&&A<=13)&&W==1) {
        if(A!=8)return {{0,true,1},{(1u<<A)-1,true,0},{(1u<<(A-1))-1,false,0},{0,false,0},{(1u<<A)-1,false,1}};
    }
    if constexpr(A==8 && W==1)
        return {{0,true,1},{255,true,0},{127,false,0},{0,false,0},
                {255,false,1},{128,true,1},{128,false,0},{127,true,0}};
    static_assert((A==2&&W==2)||((A>=8&&A<=13)&&W==1));
    return {{0,true,1},{1,true,2},{2,true,3},{3,true,0},
            {0,false,2},{1,false,3},{2,false,0},{3,false,1},
            {0,true,2},{0,true,1},{1,true,1},{2,true,2},
            {0,false,0},{1,false,0},{2,false,0},{3,false,0}};
}
template<std::size_t A,std::size_t W>
void saveTrace(const std::filesystem::path &root,const TFHEpp::Key<L> &secret,
               const std::string &filename)
{
    if(std::filesystem::exists(root/filename))return;
    auto fixture=std::make_unique<Fixture<A,W>>();
    for(std::size_t w=0;w<W;++w)for(std::size_t i=0;i<(std::size_t{1}<<A);++i){
        TFHEpp::Polynomial<L> message{};message[0]=(i>>w)&1?L::μ:-L::μ;
        TFHEpp::trlweSymEncrypt<L>(fixture->initial[w][i],message,secret);
    }
    for(const auto request:trace<A,W>()){
        typename Fixture<A,W>::RAM::Address address;
        typename Fixture<A,W>::RAM::Word word;
        TFHEpp::TLWE<L> enable;
        for(std::size_t i=0;i<A;++i)
            TFHEpp::tlweSymEncrypt<L>(address[i],(request.address>>i)&1?L::μ:-L::μ,secret);
        for(std::size_t i=0;i<W;++i)
            TFHEpp::tlweSymEncrypt<L>(word[i],(request.data>>i)&1?L::μ:-L::μ,secret);
        TFHEpp::tlweSymEncrypt<L>(enable,request.enable?L::μ:-L::μ,secret);
        fixture->addresses.push_back(address);fixture->enables.push_back(enable);fixture->data.push_back(word);
    }
    save(root/filename,*fixture);
}
} // namespace ckks_ram
