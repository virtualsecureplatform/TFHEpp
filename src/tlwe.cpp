#include <randen.h>

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
TLWE<P> tlweSymEncrypt(const typename P::T p, const double α,
                       const std::array<typename P::T, P::n> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    TLWE<P> res = {};
    res[P::n] = ModularGaussian<P>(p, α);
    for (int i = 0; i < P::n; i++) {
        res[i] = Torusdist(generator);
        res[P::n] += res[i] * key[i];
    }
    return res;
}
#define INST(P)                                \
    template TLWE<P> tlweSymEncrypt<P>(        \
        const typename P::T p, const double α, \
        const std::array<typename P::T, P::n> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::n];
    for (int i = 0; i < P::n; i++) phase -= c[i] * key[i];
    bool res =
        static_cast<typename std::make_signed<typename P::T>::type>(phase) > 0;
    return res;
}
#define INST(P) \
    template bool tlweSymDecrypt<P>(const TLWE<P> &c, const Key<P> &key)
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

vector<TLWE<lvl1param>> bootsSymEncryptHalf(const vector<uint8_t> &p,
                                            const SecretKey &sk)
{
    vector<TLWE<lvl1param>> c(p.size());
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncrypt<lvl1param>(
            p[i] ? lvl1param::μ / 2 : -(lvl1param::μ / 2), lvl1param::α,
            sk.key.lvl1);
    return c;
}

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
