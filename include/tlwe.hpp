#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "key.hpp"
#include "params.hpp"
#include "utils.hpp"

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

template <class P>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::k * P::n];
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++)
            phase -= c[k * P::n + i] * key[k * P::n + i];
    typename P::T res = static_cast<typename P::T>(std::round(phase / P::Δ)) %
                        (2 * P::plain_modulus);
    return res;
}

template <class P>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const Key<P> &key)
{
    vector<TLWE<P>> c(p.size());
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncrypt<P>(p[i] ? P::μ : -P::μ, P::α, key);
    return c;
}

template <class P = lvl1param>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const SecretKey &sk)
{
    return bootsSymEncrypt<P>(p, sk.key.get<P>());
}

template <class P>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const Key<P> &key)
{
    vector<uint8_t> p(c.size());
    for (int i = 0; i < c.size(); i++) p[i] = tlweSymDecrypt<P>(c[i], key);
    return p;
}

template <class P = lvl1param>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const SecretKey &sk)
{
    return bootsSymDecrypt<P>(c, sk.key.get<P>());
}

#include "externs/tlwe.hpp"
}  // namespace TFHEpp