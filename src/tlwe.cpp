#ifdef USE_RANDEN
#import <randen.h>
#endif

#import <array>
#import <cstdint>
#import <key.hpp>
#import <limits>
#import <params.hpp>
#import <random>
#import <tlwe.hpp>
#import <type_traits>
#import <vector>

namespace TFHEpp {

#define INST(P) \
    template TLWE<P> tlweSymEncrypt<P>(const typename P::T p, const Key<P> &key)
#undef INST

#define INST(P)                                                  \
    template TLWE<P> tlweSymIntEncrypt<P>(const typename P::T p, \
                                          const Key<P> &key)
#undef INST

#define INST(P) \
    template bool tlweSymDecrypt<P>(const TLWE<P> &c, const Key<P> &key)
#undef INST

#define INST(P)                                                   \
    template typename P::T tlweSymIntDecrypt<P>(const TLWE<P> &c, \
                                                const Key<P> &key)
#undef INST

#define INST(P)                                       \
    template std::vector<TLWE<P>> bootsSymEncrypt<P>( \
        const std::vector<uint8_t> &p, const SecretKey &sk)
#undef INST

#define INST(P)                                       \
    template std::vector<uint8_t> bootsSymDecrypt<P>( \
        const std::vector<TLWE<P>> &c, const SecretKey &sk)
#undef INST
}  // namespace TFHEpp
