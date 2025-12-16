#pragma once
#include"../tlwe.hpp"

namespace TFHEpp{
#define INST(P)                                            \
    extern template void tlweSymEncrypt<P>(TLWE<P> & res,  \
                                           const typename P::T p, \
                                           const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                               \
    extern template void tlweSymIntEncrypt<P>(TLWE<P> & res, \
                                              const typename P::T p, \
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

#define INST(P)                                                       \
    extern template void bootsSymEncrypt<P>(std::vector<TLWE<P>> & c,  \
                                            const std::vector<uint8_t> &p, \
                                            const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                       \
    extern template std::vector<uint8_t> bootsSymDecrypt<P>( \
        const std::vector<TLWE<P>> &c, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}