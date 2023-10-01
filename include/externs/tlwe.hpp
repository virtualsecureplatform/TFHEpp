#pragma once
#include"../tlwe.hpp"

namespace TFHEpp{
#define INST(P)                                                               \
    extern template TLWE<P> tlweSymEncrypt<P>(const typename P::T p, \
                                       const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                                  \
    extern template TLWE<P> tlweSymIntEncrypt<P>(const typename P::T p, \
                                          const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P) \
    extern template bool tlweSymDecrypt<P>(const TLWE<P> &c, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                                   \
    extern template typename P::T tlweSymIntDecrypt<P>(const TLWE<P> &c, \
                                                const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                       \
    extern template std::vector<TLWE<P>> bootsSymEncrypt<P>( \
        const std::vector<uint8_t> &p, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                       \
    extern template std::vector<uint8_t> bootsSymDecrypt<P>( \
        const std::vector<TLWE<P>> &c, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}