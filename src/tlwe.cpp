#ifdef USE_RANDEN
#include <randen.h>
#endif

#include <array>
#include <cstdint>
#include <key.hpp>
#include <limits>
#include <params.hpp>
#include <random>
#include <tlwe.hpp>
#include <type_traits>
#include <vector>

namespace TFHEpp {

#define INST(P) \
    template TLWE<P> tlweSymEncrypt<P>(const typename P::T p, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                                  \
    template TLWE<P> tlweSymIntEncrypt<P>(const typename P::T p, \
                                          const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P) \
    template bool tlweSymDecrypt<P>(const TLWE<P> &c, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                                   \
    template typename P::T tlweSymIntDecrypt<P>(const TLWE<P> &c, \
                                                const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                       \
    template std::vector<TLWE<P>> bootsSymEncrypt<P>( \
        const std::vector<uint8_t> &p, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

#define INST(P)                                       \
    template std::vector<uint8_t> bootsSymDecrypt<P>( \
        const std::vector<TLWE<P>> &c, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}  // namespace TFHEpp
