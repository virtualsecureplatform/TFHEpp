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
TLWE<P> tlweSymEncrypt(const typename P::T p, const uint η, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(0, P::q - 1);
    TLWE<P> res = {};
    res[P::k * P::n] =
        p + CenteredBinomial<P>(η)
        << (std::numeric_limits<typename P::T>::digits - P::qbit);
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++) {
            res[k * P::n + i] =
                Torusdist(generator)
                << (std::numeric_limits<typename P::T>::digits - P::qbit);
            res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
        }
    return res;
}

template <class P>
TLWE<P> tlweSymEncrypt(const typename P::T p, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return tlweSymEncrypt<P>(p, P::α, key);
    else
        return tlweSymEncrypt<P>(p, P::η, key);
}

template <class P>
TLWE<P> tlweSymEncrypt(const typename P::T p, const SecretKey &sk)
{
    return tlweSymEncrypt<P>(p, sk.key.get<P>());
}

template <class P, uint plain_modulus = P::plain_modulus>
TLWE<P> tlweSymIntEncrypt(const typename P::T p, const double α,
                          const Key<P> &key)
{
    constexpr double Δ =
        std::pow(2.0, std::numeric_limits<typename P::T>::digits) /
        plain_modulus;
    return tlweSymEncrypt<P>(static_cast<typename P::T>(p * Δ), α, key);
}

template <class P, uint plain_modulus = P::plain_modulus>
TLWE<P> tlweSymIntEncrypt(const typename P::T p, const uint η,
                          const Key<P> &key)
{
    constexpr double Δ =
        std::pow(2.0, std::numeric_limits<typename P::T>::digits) /
        plain_modulus;
    return tlweSymEncrypt<P>(static_cast<typename P::T>(p * Δ), η, key);
}

template <class P, uint plain_modulus = P::plain_modulus>
TLWE<P> tlweSymIntEncrypt(const typename P::T p, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return tlweSymIntEncrypt<P, plain_modulus>(p, P::α, key);
    else
        return tlweSymIntEncrypt<P, plain_modulus>(p, P::η, key);
}

template <class P, uint plain_modulus = P::plain_modulus>
TLWE<P> tlweSymIntEncrypt(const typename P::T p, const SecretKey &sk)
{
    return tlweSymIntEncrypt<P, plain_modulus>(p, sk.key.get<P>());
}

template <class P>
typename P::T tlweSymPhase(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::k * P::n];
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++)
            phase -= c[k * P::n + i] * key[k * P::n + i];
    return phase;
}

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = tlweSymPhase<P>(c, key);
    bool res =
        static_cast<typename std::make_signed<typename P::T>::type>(phase) > 0;
    return res;
}

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const SecretKey &sk)
{
    return tlweSymDecrypt<P>(c, sk.key.get<P>());
}

template <class P, uint plain_modulus = P::plain_modulus>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    constexpr double Δ =
        2 *
        static_cast<double>(
            1ULL << (std::numeric_limits<typename P::T>::digits - 1)) /
        plain_modulus;
    const typename P::T phase = tlweSymPhase<P>(c, key);
    typename P::T res = static_cast<typename P::T>(std::round(phase / Δ));
    return res >= plain_modulus / 2 ? res - plain_modulus : res;
}

template <class P, uint plain_modulus = P::plain_modulus>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const SecretKey &sk)
{
    return tlweSymIntDecrypt<P, plain_modulus>(c, sk.key.get<P>());
}

template <class P, std::make_signed_t<typename P::T> μ>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const Key<P> &key)
{
    vector<TLWE<P>> c(p.size());
#pragma omp parallel for
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncrypt<P>(p[i] ? μ : -μ, key);
    return c;
}

template <class P>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const Key<P> &key)
{
    return bootsSymEncrypt<P, P::μ>(p, key);
}

template <class P = lvl1param>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const SecretKey &sk)
{
    return bootsSymEncrypt<P>(p, sk.key.get<P>());
}

template <class P = lvl1param, std::make_signed_t<typename P::T> μ>
std::vector<TLWE<P>> bootsSymEncrypt(const std::vector<uint8_t> &p,
                                     const SecretKey &sk)
{
    return bootsSymEncrypt<P, μ>(p, sk.key.get<P>());
}

template <class P>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const Key<P> &key)
{
    vector<uint8_t> p(c.size());
#pragma omp parallel for
    for (int i = 0; i < c.size(); i++) p[i] = tlweSymDecrypt<P>(c[i], key);
    return p;
}

template <class P = lvl1param>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const SecretKey &sk)
{
    return bootsSymDecrypt<P>(c, sk.key.get<P>());
}

template <class P>
void TLWEAdd(TLWE<P> &res, const TLWE<P> &a, const TLWE<P> &b)
{
    for (int i = 0; i <= P::k * P::n; i++) res[i] = a[i] + b[i];
}

}  // namespace TFHEpp