#pragma once

#include "key.hpp"
#include "mulfft.hpp"
#include "params.hpp"

namespace TFHEpp {
template <class P>
void trlweSymEncryptZero(TRLWE<P> &c, const double α, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    for (typename P::T &i : c[P::k]) i = ModularGaussian<P>(0, α);
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
}

template <class P>
void trlweSymEncryptZero(TRLWE<P> &c, const uint η, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
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
}

template <class P>
void trlweSymEncryptZero(TRLWE<P> &c, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trlweSymEncryptZero<P>(c, P::α, key);
    else
        trlweSymEncryptZero<P>(c, P::η, key);
}

template <class P>
void trlweSymEncryptZero(TRLWERAINTT<P> &c, const uint η, const Key<P> &key)
{
    static_assert(P::q == raintt::P);
    static_assert(P::qbit == raintt::wordbits);
    std::uniform_int_distribution<typename P::T> Torusdist(0, P::q - 1);
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    c = {};
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
}

template <class P>
void trlweSymEncrypt(TRLWE<P> &c, const std::array<typename P::T, P::n> &p,
                     const double α, const Key<P> &key)
{
    trlweSymEncryptZero<P>(c, α, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
}

template <class P>
void trlweSymEncrypt(TRLWE<P> &c, const std::array<typename P::T, P::n> &p,
                     const uint η, const Key<P> &key)
{
    trlweSymEncryptZero<P>(c, η, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
}

template <class P>
void trlweSymEncrypt(TRLWE<P> &c, const std::array<typename P::T, P::n> &p,
                     const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trlweSymEncrypt<P>(c, p, P::α, key);
    else
        trlweSymEncrypt<P>(c, p, P::η, key);
}

template <class P, bool modswitch = false>
void trlweSymEncrypt(TRLWERAINTT<P> &c, const Polynomial<P> &p, const uint η,
                     const Key<P> &key)
{
    trlweSymEncryptZero<P>(c, η, key);
    PolynomialRAINTT<P> pntt;
    raintt::TwistINTT<typename P::T, P::nbit, modswitch>(
        pntt, p, (*raintttable)[1], (*raintttwist)[1]);
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    for (int i = 0; i < P::n; i++)
        if ((i & ((1 << remainder) - 1)) > 1)
            pntt[i] = raintt::MulSREDC(pntt[i], raintt::R2);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] = raintt::AddMod(pntt[i], c[P::k][i]);
}

template <class P>
void trlweSymIntEncrypt(TRLWE<P> &c, const std::array<typename P::T, P::n> &p,
                        const double α, const Key<P> &key)
{
    trlweSymEncryptZero<P>(c, α, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
}

template <class P>
void trlweSymIntEncrypt(TRLWE<P> &c, const std::array<typename P::T, P::n> &p,
                        const uint η, const Key<P> &key)
{
    trlweSymEncryptZero<P>(c, η, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
}

template <class P>
void trlweSymIntEncrypt(TRLWE<P> &c, const std::array<typename P::T, P::n> &p,
                        const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trlweSymIntEncrypt<P>(c, p, P::α, key);
    else
        trlweSymIntEncrypt<P>(c, p, P::η, key);
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

template <class P, uint plain_modulus = P::plain_modulus>
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

    const double Δ = std::pow(2.0, std::numeric_limits<typename P::T>::digits) /
                     plain_modulus;
    Polynomial<P> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(std::round(phase[i] / Δ)) %
               plain_modulus;
    return p;
}

template <class P, uint plain_modulus = P::plain_modulus>
Polynomial<P> trlweSymIntDecrypt(const TRLWE<P> &c, const SecretKey &sk)
{
    return trlweSymIntDecrypt<P, plain_modulus>(c, sk.key.get<P>());
}

template <class P, class... Args>
void TRLWEAdd(TRLWE<P> &res, const TRLWE<P> &first, const Args &...rest)
{
    for (int i = 0; i < (P::k + 1) * P::n; i++) {
        // A binary fold expression sums all corresponding elements at once.
        res[0][i] = (first[0][i] + ... + rest[0][i]);
    }
}

template <class P, class... Args>
void TRLWESub(TRLWE<P> &res, const TRLWE<P> &first, const Args &...rest)
{
    for (int i = 0; i < (P::k + 1) * P::n; i++) {
        // A binary fold expression subtracts all corresponding elements at
        // once.
        res[0][i] = (first[0][i] - ... - rest[0][i]);
    }
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
