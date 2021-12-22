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
#include <utils.hpp>
#include <vector>

namespace TFHEpp {

template <class P>
TLWE<P> tlweSymEncrypt(const typename P::T p, const double α, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    TLWE<P> res = {};
    res[P::k * P::n] = ModularGaussian<P>(p, α);
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++) {
            res[k * P::n + i] = Torusdist(generator);
            res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
        }
    return res;
}
#define INST(P)                                                               \
    template TLWE<P> tlweSymEncrypt<P>(const typename P::T p, const double α, \
                                       const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
TLWE<P> tlweSymIntEncrypt(const typename P::T p, const double α,
                          const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    TLWE<P> res = {};
    res[P::k * P::n] =
        ModularGaussian<P>(static_cast<typename P::T>(p * P::Δ), α);
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++) {
            res[k * P::n + i] = Torusdist(generator);
            res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
        }
    return res;
}
#define INST(P)                                                  \
    template TLWE<P> tlweSymIntEncrypt<P>(const typename P::T p, \
                                          const double α, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::k * P::n];
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++)
            phase -= c[k * P::n + i] * key[k * P::n + i];
    bool res =
        static_cast<typename std::make_signed<typename P::T>::type>(phase) > 0;
    return res;
}
#define INST(P) \
    template bool tlweSymDecrypt<P>(const TLWE<P> &c, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::k * P::n];
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++)
            phase -= c[k * P::n + i] * key[k * P::n + i];
    typename P::T res =
        static_cast<typename P::T>(std::round(phase / P::Δ)) % P::plain_modulus;
    return res;
}
#define INST(P)                                                   \
    template typename P::T tlweSymIntDecrypt<P>(const TLWE<P> &c, \
                                                const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const SecretKey &sk)
{
    vector<TLWE<P>> c(p.size());
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncrypt<P>(p[i] ? P::μ : -P::μ, P::α, sk.key.get<P>());
    return c;
}
#define INST(P)                                       \
    template std::vector<TLWE<P>> bootsSymEncrypt<P>( \
        const std::vector<uint8_t> &p, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const SecretKey &sk)
{
    vector<uint8_t> p(c.size());
    for (int i = 0; i < c.size(); i++)
        p[i] = tlweSymDecrypt<P>(c[i], sk.key.get<P>());
    return p;
}
#define INST(P)                                       \
    template std::vector<uint8_t> bootsSymDecrypt<P>( \
        const std::vector<TLWE<P>> &c, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}  // namespace TFHEpp
