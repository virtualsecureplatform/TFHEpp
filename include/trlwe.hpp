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
        0, std::numeric_limits<typename P::T>::max());
    alignas(64) TRLWE<P> c;
    for (typename P::T &i : c[P::k])
        i = (CenteredBinomial<P>(η) << std::numeric_limits<P>::digits) / P::q;
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        alignas(64) std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        alignas(64) Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
    return c;
}

template <class P>
TRLWE<P> trlweSymEncryptZero(const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncryptZero<P>(P::α, key);
    else
        return trlweSymEncryptZero<P>(P::η, key);
}

template <class P>
TRLWERAINTT<P> trlwerainttSymEncryptZero(const uint η, const Key<P> &key)
{
    static_assert(P::q == raintt::P);
    static_assert(P::qbit == raintt::wordbits);
    std::uniform_int_distribution<typename P::T> Torusdist(0, P::q - 1);
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    TRLWERAINTT<P> c = {};
    {
        Polynomial<P> b;
        for (typename P::T &i : b) i = CenteredBinomial<P>(η);
        raintt::TwistINTT<typename P::T, P::nbit, false>(
            c[P::k], b, (*raintttable)[1], (*raintttwist)[1]);
        for (int i = 0; i < P::n; i++)
            if ((i & ((1 << remainder) - 1)) > 1)
                c[P::k][i] = raintt::MulSREDC(c[P::k][i], raintt::R2);
    }
    for (int k = 0; k < P::k; k++) {
        for (typename raintt::DoubleSWord &i : c[k]) i = Torusdist(generator);
        PolynomialRAINTT<P> partkeyraintt;
        {
            Polynomial<P> partkey;
            for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
            raintt::TwistINTT<typename P::T, P::nbit, false>(
                partkeyraintt, partkey, (*raintttable)[1], (*raintttwist)[1]);
            for (int i = 0; i < P::n; i++)
                if ((i & ((1 << remainder) - 1)) > 1)
                    partkeyraintt[i] =
                        raintt::MulSREDC(partkeyraintt[i], raintt::R3);
                else
                    partkeyraintt[i] =
                        raintt::MulSREDC(partkeyraintt[i], raintt::R2);
        }
        for (int i = 0; i < P::n; i++)
            c[P::k][i] = raintt::AddMod(
                c[P::k][i], raintt::MulSREDC(c[k][i], partkeyraintt[i]));
    }
    return c;
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
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p, const uint η,
                         const Key<P> &key)
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
        return trlweSymEncrypt<P>(p, P::α, key);
    else
        return trlweSymEncrypt<P>(p, P::η, key);
}

template <class P, bool modswitch = false>
TRLWERAINTT<P> trlwerainttSymEncrypt(const Polynomial<P> &p, const uint η,
                                     const Key<P> &key)
{
    TRLWERAINTT<P> c = trlwerainttSymEncryptZero<P>(η, key);
    PolynomialRAINTT<P> pntt;
    raintt::TwistINTT<typename P::T, P::nbit, modswitch>(
        pntt, p, (*raintttable)[1], (*raintttwist)[1]);
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    for (int i = 0; i < P::n; i++)
        if ((i & ((1 << remainder) - 1)) > 1)
            pntt[i] = raintt::MulSREDC(pntt[i], raintt::R2);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] = raintt::AddMod(pntt[i], c[P::k][i]);
    return c;
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
        return trlweSymIntEncrypt<P>(p, P::α, key);
    else
        return trlweSymIntEncrypt<P>(p, P::η, key);
}

template <class P>
Polynomial<P> trlwePhase(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }
    return phase;
}

template <class P>
std::array<bool, P::n> trlweSymDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(c, key);

    std::array<bool, P::n> p;
    if constexpr (hasq<P>) {
        for (int i = 0; i < P::n; i++) p[i] = (phase[i] % P::q) < P::q / 2;
    }
    else
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