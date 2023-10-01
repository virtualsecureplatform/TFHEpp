#pragma once

#include "mulfft.hpp"
#include "params.hpp"

namespace TFHEpp {
template <class P>
TRLWE<P> trlweSymEncryptZero(const double α, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    TRLWE<P> c;
    for (typename P::T &i : c[P::k]) i = ModularGaussian<P>(0, α);
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
    return c;
}

template <class P>
TRLWE<P> trlweSymEncryptZero(const uint η, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, P::q-1);
    TRLWE<P> c;
    for (typename P::T &i : c[P::k]) i = CenteredBinomial<P>(η)<<(std::numeric_limits<typename P::T>::digits-P::qbit);
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator)<<(std::numeric_limits<typename P::T>::digits-P::qbit);
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
    return c;
}

template <class P>
TRLWE<P> trlweSymEncryptZero(const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncryptZero<P>(P::α,key); 
    else
        return trlweSymEncryptZero<P>(P::η,key); 
}

template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const double α, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}

template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const uint η, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(η, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}

template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncrypt<P>(p, P::α,key); 
    else
        return trlweSymEncrypt<P>(p, P::η,key); 
}

template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const double α, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
    return c;
}

template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const uint η, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(η, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
    return c;
}

template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymIntEncrypt<P>(p, P::α,key); 
    else
        return trlweSymIntEncrypt<P>(p, P::η,key); 
}

template <class P>
std::array<bool, P::n> trlweSymDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }

    std::array<bool, P::n> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename std::make_signed<typename P::T>::type>(
                   phase[i]) > 0;
    return p;
}

template <class P>
Polynomial<P> trlweSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }

    Polynomial<P> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(std::round(phase[i] / P::Δ)) %
               P::plain_modulus;
    return p;
}

template <class P>
void SampleExtractIndex(TLWE<P> &tlwe, const TRLWE<P> &trlwe, const int index)
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            tlwe[k * P::n + i] = trlwe[k][index - i];
        for (int i = index + 1; i < P::n; i++)
            tlwe[k * P::n + i] = -trlwe[k][P::n + index - i];
    }
    tlwe[P::k * P::n] = trlwe[P::k][index];
}

template <class P>
void InvSampleExtractIndex(TRLWE<P> &trlwe, const TLWE<P> &tlwe,
                           const int index)
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            trlwe[k][index - i] = tlwe[k * P::n + i];
        for (int i = index + 1; i < P::n; i++)
            trlwe[k][P::n + index - i] = -tlwe[k * P::n + i];
    }
    trlwe[P::k] = {};
    trlwe[P::k][index] = tlwe[P::k * P::n];
}

}  // namespace TFHEpp