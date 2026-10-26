#pragma once
#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include "params.hpp"

namespace TFHEpp {
template<class P> inline double TFHEParameterNoiseValue()
{
    if constexpr(requires{P::α;})return P::α;
    else if constexpr(requires{P::η;})return double(P::η);
    else return TFHEParameterNoiseValue<typename P::targetP>();
}
template<class P> inline auto TFHELevelProfile()
{
    std::uint64_t noise;
    if constexpr(P::errordist==ErrorDistribution::ModularGaussian)
        noise=std::bit_cast<std::uint64_t>(double(P::α));
    else noise=P::η;
    std::array<std::uint64_t,14> result{P::n,P::k,std::numeric_limits<typename P::T>::digits,
        std::uint64_t(P::errordist),noise,0,0,0,0,0,std::numeric_limits<typename P::T>::digits,
        std::uint64_t(P::key_value_min),std::uint64_t(P::key_value_max),0};
    if constexpr(requires{P::l;})result[5]=P::l;
    if constexpr(requires{P::Bgbit;})result[6]=P::Bgbit;
    if constexpr(requires{P::lₐ;})result[7]=P::lₐ;
    if constexpr(requires{P::Bgₐbit;})result[8]=P::Bgₐbit;
    if constexpr(requires{P::q;})result[9]=std::uint64_t(P::q);
    if constexpr(requires{P::qbit;})result[10]=P::qbit;
    if constexpr(requires{P::ell;})result[13]=P::ell;
    return result;
}
inline auto TFHEArchiveNoiseProfile()
{
    std::array<std::uint64_t,54> profile{};
    std::size_t next=0;
    for(auto value:TFHELevelProfile<lvl0param>())profile[next++]=value;
    for(auto value:TFHELevelProfile<lvl1param>())profile[next++]=value;
    for(auto value:TFHELevelProfile<lvl2param>())profile[next++]=value;
    for(auto value:std::array<std::uint64_t,9>{lvl10param::t,lvl10param::basebit,std::bit_cast<std::uint64_t>(TFHEParameterNoiseValue<lvl10param>()),
        lvl20param::t,lvl20param::basebit,std::bit_cast<std::uint64_t>(TFHEParameterNoiseValue<lvl20param>()),
        lvl21param::t,lvl21param::basebit,std::bit_cast<std::uint64_t>(TFHEParameterNoiseValue<lvl21param>())})profile[next++]=value;
#ifdef USE_BLAKE3
    profile[next++]=1;
#elif defined(USE_RANDEN)
    profile[next++]=2;
#else
    profile[next++]=3;
#endif
    profile[next++]=1;
    profile[next++]=2;
    return profile;
}
} // namespace TFHEpp
