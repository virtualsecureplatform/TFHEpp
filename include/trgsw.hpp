#pragma once

#include <array>
#include <cstdint>

#include "mulfft.hpp"
#include "params.hpp"
#include "trlwe.hpp"

namespace TFHEpp {

template <class P>
constexpr typename P::T offsetgen()
{
    typename P::T offset = 0;
    for (int i = 1; i <= P::l; i++)
        offset += P::Bg / 2 *
                  (static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                            i * P::Bgbit));
    return offset;
}

#ifdef USE_OPTIMAL_DECOMPOSITION
template <class P>
inline void Decomposition(DecomposedPolynomial<P> &decpoly,
                          const Polynomial<P> &poly, typename P::T randbits)
{
    // https://eprint.iacr.org/2021/1161
    constexpr typename P::T roundoffset =
        static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits - P::l * P::Bgbit -
                 1);
    constexpr typename P::T mask =
        static_cast<typename P::T>((static_cast<typename P::T>(1) << P::Bgbit) - 1);
    constexpr typename P::T Bgl = static_cast<typename P::T>(1) << (P::l * P::Bgbit);

    Polynomial<P> K;
    for (int i = 0; i < P::n; i++) {
        K[i] = (poly[i] + roundoffset) >>
               (std::numeric_limits<typename P::T>::digits - P::l * P::Bgbit);
        if (K[i] > Bgl / 2)
            K[i] -= Bgl;
        else if (K[i] == Bgl / 2) {
            if (randbits & 1) K[i] -= Bgl;
            randbits >>= 1;
        }
    }
    for (int l = 0; l < P::l; l++) {
        for (int i = 0; i < P::n; i++) {
            // https://github.com/zama-ai/tfhe-rs/blob/06b700f9042eb5dfbaf073bb6b7f71bff4be1c2f/tfhe/src/core_crypto/commons/math/decomposition/iter.rs#L117-L124
            auto &ki = decpoly[P::l - l - 1][i];
            ki = K[i] & mask;
            K[i] >>= P::Bgbit;
            // if((ki>P::Bg/2)||((ki==P::Bg/2)&&((K[i]&mask)>=P::Bg/2))){
            //     ki -= P::Bg;
            //     K[i] += 1;
            // }
            const auto carry = (((ki - 1) | K[i]) & ki) >> (P::Bgbit - 1);
            K[i] += carry;
            ki -= carry << P::Bgbit;
        }
    }
}
#endif

template <class P>
inline void Decomposition(DecomposedPolynomial<P> &decpoly,
                          const Polynomial<P> &poly)
{
    constexpr typename P::T offset = offsetgen<P>();
    constexpr typename P::T roundoffset =
        static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits - P::l * P::Bgbit -
                 1);
    constexpr typename P::T mask =
        static_cast<typename P::T>((static_cast<typename P::T>(1) << P::Bgbit) - 1);
    constexpr typename P::T halfBg = (static_cast<typename P::T>(1) << (P::Bgbit - 1));

    for (int i = 0; i < P::n; i++) {
        for (int l = 0; l < P::l; l++)
            decpoly[l][i] = (((poly[i] + offset + roundoffset) >>
                              (std::numeric_limits<typename P::T>::digits -
                               (l + 1) * P::Bgbit)) &
                             mask) -
                            halfBg;
    }
}

template <class P>
constexpr typename P::T nonceoffsetgen()
{
    typename P::T offset = 0;
    for (int i = 1; i <= P::lₐ; i++)
        offset += P::Bgₐ / 2 *
                  (static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                            i * P::Bgₐbit));
    return offset;
}

template <class P>
inline void NonceDecomposition(DecomposedNoncePolynomial<P> &decpoly,
                               const Polynomial<P> &poly)
{
    constexpr typename P::T offset = nonceoffsetgen<P>();
    constexpr typename P::T roundoffset =
        static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                 P::lₐ * P::Bgₐbit - 1);
    constexpr typename P::T mask =
        static_cast<typename P::T>((static_cast<typename P::T>(1) << P::Bgₐbit) - 1);
    constexpr typename P::T halfBg = (static_cast<typename P::T>(1) << (P::Bgₐbit - 1));

    for (int i = 0; i < P::n; i++) {
        for (int l = 0; l < P::lₐ; l++)
            decpoly[l][i] = (((poly[i] + offset + roundoffset) >>
                              (std::numeric_limits<typename P::T>::digits -
                               (l + 1) * P::Bgₐbit)) &
                             mask) -
                            halfBg;
    }
}

// Double Decomposition (bivariate representation) for external product
// Decomposes each coefficient a into l*l̅ components such that:
// a ≈ Σᵢ Σⱼ aᵢⱼ * Bg^(l-i) * B̅g^(l̅-j)
// When l̅=1 (j=0 only), this reduces to standard decomposition.
template <class P>
constexpr typename P::T ddoffsetgen()
{
    typename P::T offset = 0;
    for (int i = 1; i <= P::l; i++)
        for (int j = 0; j < P::l̅; j++)
            offset += (static_cast<typename P::T>(P::Bg) / 2) *
                      (static_cast<typename P::T>(1)
                       << (std::numeric_limits<typename P::T>::digits -
                           i * P::Bgbit - j * P::B̅gbit));
    return offset;
}

template <class P>
inline void DoubleDecomposition(DecomposedPolynomialDD<P> &decpoly,
                                const Polynomial<P> &poly)
{
    constexpr typename P::T offset = ddoffsetgen<P>();
    // Remaining bits after decomposition
    constexpr int remaining_bits = std::numeric_limits<typename P::T>::digits -
                                   P::l * P::Bgbit - P::l̅ * P::B̅gbit;
    // roundoffset is 0 if no remaining bits, otherwise 2^(remaining_bits-1)
    constexpr typename P::T roundoffset =
        remaining_bits > 0
            ? (static_cast<typename P::T>(1) << (remaining_bits - 1))
            : static_cast<typename P::T>(0);
    constexpr typename P::T maskBg =
        static_cast<typename P::T>((static_cast<typename P::T>(1) << P::Bgbit) - 1);
    constexpr typename P::T halfBg = (static_cast<typename P::T>(1) << (P::Bgbit - 1));

    for (int n = 0; n < P::n; n++) {
        typename P::T a = poly[n] + offset + roundoffset;
        for (int i = 0; i < P::l; i++) {
            for (int j = 0; j < P::l̅; j++) {
                // Shift to get the (i,j)-th digit in base Bg (after B̅g grouping)
                // When l̅=1 (j=0 only), this reduces to standard decomposition
                const int shift = std::numeric_limits<typename P::T>::digits -
                                  (i + 1) * P::Bgbit - j * P::B̅gbit;
                auto decomp_val =
                    static_cast<std::make_signed_t<typename P::T>>(
                        ((a >> shift) & maskBg) - halfBg);
                // For 128-bit types, shift left by 64 so TwistIFFT (which uses
                // top 64 bits) gets the correct small integer value
                if constexpr (std::is_same_v<typename P::T, __uint128_t>)
                    decpoly[i * P::l̅ + j][n] =
                        static_cast<typename P::T>(decomp_val) << 64;
                else
                    decpoly[i * P::l̅ + j][n] = decomp_val;
            }
        }
    }
}

template <class P>
constexpr typename P::T nonceddoffsetgen()
{
    typename P::T offset = 0;
    for (int i = 1; i <= P::lₐ; i++)
        for (int j = 0; j < P::l̅ₐ; j++)
            offset += (static_cast<typename P::T>(P::Bgₐ) / 2) *
                      (static_cast<typename P::T>(1)
                       << (std::numeric_limits<typename P::T>::digits -
                           i * P::Bgₐbit - j * P::B̅gₐbit));
    return offset;
}

template <class P>
inline void NonceDoubleDecomposition(DecomposedNoncePolynomialDD<P> &decpoly,
                                     const Polynomial<P> &poly)
{
    constexpr typename P::T offset = nonceddoffsetgen<P>();
    // Remaining bits after decomposition
    constexpr int remaining_bits = std::numeric_limits<typename P::T>::digits -
                                   P::lₐ * P::Bgₐbit - P::l̅ₐ * P::B̅gₐbit;
    // roundoffset is 0 if no remaining bits, otherwise 2^(remaining_bits-1)
    constexpr typename P::T roundoffset =
        remaining_bits > 0
            ? (static_cast<typename P::T>(1) << (remaining_bits - 1))
            : static_cast<typename P::T>(0);
    constexpr typename P::T maskBg =
        static_cast<typename P::T>((static_cast<typename P::T>(1) << P::Bgₐbit) - 1);
    constexpr typename P::T halfBg = (static_cast<typename P::T>(1) << (P::Bgₐbit - 1));

    for (int n = 0; n < P::n; n++) {
        typename P::T a = poly[n] + offset + roundoffset;
        for (int i = 0; i < P::lₐ; i++) {
            for (int j = 0; j < P::l̅ₐ; j++) {
                // Shift to get the (i,j)-th digit
                // When l̅ₐ=1 (j=0 only), this reduces to standard decomposition
                const int shift = std::numeric_limits<typename P::T>::digits -
                                  (i + 1) * P::Bgₐbit - j * P::B̅gₐbit;
                auto decomp_val =
                    static_cast<std::make_signed_t<typename P::T>>(
                        ((a >> shift) & maskBg) - halfBg);
                // For 128-bit types, shift left by 64 so TwistIFFT (which uses
                // top 64 bits) gets the correct small integer value
                if constexpr (std::is_same_v<typename P::T, __uint128_t>)
                    decpoly[i * P::l̅ₐ + j][n] =
                        static_cast<typename P::T>(decomp_val) << 64;
                else
                    decpoly[i * P::l̅ₐ + j][n] = decomp_val;
            }
        }
    }
}

template <class P>
void Decomposition(DecomposedPolynomialNTT<P> &decpolyntt,
                   const Polynomial<P> &poly)
{
    DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, poly);
    for (int i = 0; i < P::l; i++) TwistINTT<P>(decpolyntt[i], decpoly[i]);
}

template <class P>
void NonceDecomposition(DecomposedNoncePolynomialNTT<P> &decpolyntt,
                        const Polynomial<P> &poly)
{
    DecomposedNoncePolynomial<P> decpoly;
    NonceDecomposition<P>(decpoly, poly);
    for (int i = 0; i < P::lₐ; i++) TwistINTT<P>(decpolyntt[i], decpoly[i]);
}

template <class P>
void Decomposition(DecomposedPolynomialRAINTT<P> &decpolyntt,
                   const Polynomial<P> &poly)
{
    DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, poly);
    for (int i = 0; i < P::l; i++)
        raintt::TwistINTT<typename P::T, P::nbit, false>(
            decpolyntt[i], decpoly[i], (*raintttable)[1], (*raintttwist)[1]);
}

template <class P>
void NonceDecomposition(DecomposedNoncePolynomialRAINTT<P> &decpolyntt,
                        const Polynomial<P> &poly)
{
    DecomposedNoncePolynomial<P> decpoly;
    NonceDecomposition<P>(decpoly, poly);
    for (int i = 0; i < P::lₐ; i++)
        raintt::TwistINTT<typename P::T, P::nbit, false>(
            decpolyntt[i], decpoly[i], (*raintttable)[1], (*raintttwist)[1]);
}

template <class P>
void ExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                     const TRGSWFFT<P> &trgswfft)
{
    alignas(64) PolynomialInFD<P> decpolyfft;
    alignas(64) TRLWEInFD<P> restrlwefft;
    {
        alignas(64) DecomposedNoncePolynomial<P> decpoly;
        NonceDecomposition<P>(decpoly, trlwe[0]);
        TwistIFFT<P>(decpolyfft, decpoly[0]);
        for (int m = 0; m < P::k + 1; m++)
            MulInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[0][m]);
        for (int i = 1; i < P::lₐ; i++) {
            TwistIFFT<P>(decpolyfft, decpoly[i]);
            for (int m = 0; m < P::k + 1; m++)
                FMAInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[i][m]);
        }
        for (int k = 1; k < P::k; k++) {
            NonceDecomposition<P>(decpoly, trlwe[k]);
            for (int i = 0; i < P::lₐ; i++) {
                TwistIFFT<P>(decpolyfft, decpoly[i]);
                for (int m = 0; m < P::k + 1; m++)
                    FMAInFD<P::n>(restrlwefft[m], decpolyfft,
                                  trgswfft[i + k * P::lₐ][m]);
            }
        }
    }
    alignas(64) DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, trlwe[P::k]);
    for (int i = 0; i < P::l; i++) {
        TwistIFFT<P>(decpolyfft, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
            FMAInFD<P::n>(restrlwefft[m], decpolyfft,
                          trgswfft[i + P::k * P::lₐ][m]);
    }
    for (int k = 0; k < P::k + 1; k++) TwistFFT<P>(res[k], restrlwefft[k]);
}

// External product with Double Decomposition (bivariate representation)
// Uses the full TRGSW structure with l*l̅ rows for "b" block and k*lₐ*l̅ₐ rows
// for "a" blocks
template <class P>
void ExternalProductDD(TRLWE<P> &res, const TRLWE<P> &trlwe,
                       const TRGSWFFT<P> &trgswfft)
{
    alignas(64) PolynomialInFD<P> decpolyfft;
    alignas(64) TRLWEInFD<P> restrlwefft;

    // Handle "a" polynomials (indices 0 to k-1 in TRLWE)
    // Uses NonceDoubleDecomposition with lₐ*l̅ₐ levels
    {
        alignas(64) DecomposedNoncePolynomialDD<P> decpoly;
        NonceDoubleDecomposition<P>(decpoly, trlwe[0]);
        TwistIFFT<P>(decpolyfft, decpoly[0]);
        for (int m = 0; m < P::k + 1; m++)
            MulInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[0][m]);
        for (int i = 1; i < P::lₐ * P::l̅ₐ; i++) {
            TwistIFFT<P>(decpolyfft, decpoly[i]);
            for (int m = 0; m < P::k + 1; m++)
                FMAInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[i][m]);
        }
        for (int k = 1; k < P::k; k++) {
            NonceDoubleDecomposition<P>(decpoly, trlwe[k]);
            for (int i = 0; i < P::lₐ * P::l̅ₐ; i++) {
                TwistIFFT<P>(decpolyfft, decpoly[i]);
                for (int m = 0; m < P::k + 1; m++)
                    FMAInFD<P::n>(restrlwefft[m], decpolyfft,
                                  trgswfft[i + k * P::lₐ * P::l̅ₐ][m]);
            }
        }
    }

    // Handle "b" polynomial (index k in TRLWE)
    // Uses DoubleDecomposition with l*l̅ levels
    alignas(64) DecomposedPolynomialDD<P> decpoly;
    DoubleDecomposition<P>(decpoly, trlwe[P::k]);
    for (int i = 0; i < P::l * P::l̅; i++) {
        TwistIFFT<P>(decpolyfft, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
            FMAInFD<P::n>(restrlwefft[m], decpolyfft,
                          trgswfft[i + P::k * P::lₐ * P::l̅ₐ][m]);
    }
    for (int k = 0; k < P::k + 1; k++) TwistFFT<P>(res[k], restrlwefft[k]);
}

template <class P>
void ExternalProduct(TRLWE<P> &res, const Polynomial<P> &poly,
                     const HalfTRGSWFFT<P> &halftrgswfft)
{
    alignas(64) DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, poly);
    alignas(64) PolynomialInFD<P> decpolyfft;
    // __builtin_prefetch(trgswfft[0].data());
    TwistIFFT<P>(decpolyfft, decpoly[0]);
    alignas(64) TRLWEInFD<P> restrlwefft;
    for (int m = 0; m < P::k + 1; m++)
        MulInFD<P::n>(restrlwefft[m], decpolyfft, halftrgswfft[0][m]);
    for (int i = 1; i < P::l; i++) {
        // __builtin_prefetch(trgswfft[i].data());
        TwistIFFT<P>(decpolyfft, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
            FMAInFD<P::n>(restrlwefft[m], decpolyfft, halftrgswfft[i][m]);
    }
    for (int k = 0; k < P::k + 1; k++) TwistFFT<P>(res[k], restrlwefft[k]);
}

template <class P>
void ExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                     const TRGSWRAINTT<P> &trgswntt)
{
    TRLWERAINTT<P> restrlwentt;
    {
        DecomposedNoncePolynomialRAINTT<P> decpolyntt;
        NonceDecomposition<P>(decpolyntt, trlwe[0]);
        for (int m = 0; m < P::k + 1; m++)
            for (int i = 0; i < P::n; i++)
                restrlwentt[m][i] =
                    raintt::MulSREDC(decpolyntt[0][i], trgswntt[0][m][i]);
        for (int i = 1; i < P::lₐ; i++) {
            for (int m = 0; m < P::k + 1; m++)
                for (int j = 0; j < P::n; j++)
                    restrlwentt[m][j] = raintt::AddMod(
                        restrlwentt[m][j],
                        raintt::MulSREDC(decpolyntt[i][j], trgswntt[i][m][j]));
        }
        for (int k = 1; k < P::k; k++) {
            NonceDecomposition<P>(decpolyntt, trlwe[k]);
            for (int i = 0; i < P::lₐ; i++) {
                for (int m = 0; m < P::k + 1; m++)
                    for (int j = 0; j < P::n; j++)
                        restrlwentt[m][j] = raintt::AddMod(
                            restrlwentt[m][j],
                            raintt::MulSREDC(decpolyntt[i][j],
                                             trgswntt[i + k * P::lₐ][m][j]));
            }
        }
    }
    DecomposedPolynomialRAINTT<P> decpolyntt;
    Decomposition<P>(decpolyntt, trlwe[P::k]);
    for (int i = 0; i < P::l; i++) {
        for (int m = 0; m < P::k + 1; m++)
            for (int j = 0; j < P::n; j++)
                restrlwentt[m][j] = raintt::AddMod(
                    restrlwentt[m][j],
                    raintt::MulSREDC(decpolyntt[i][j],
                                     trgswntt[i + P::k * P::lₐ][m][j]));
    }
    // if constexpr(hasq<P>) for (int k = 0; k < P::k + 1; k++)
    // raintt::TwistNTT<typename P::T,P::nbit, P::q == (1ULL<<P::qbit)>(res[k],
    // restrlwentt[k],(*raintttable)[0],(*raintttwist)[0]);
    for (int k = 0; k < P::k + 1; k++)
        raintt::TwistNTT<typename P::T, P::nbit, true>(
            res[k], restrlwentt[k], (*raintttable)[0], (*raintttwist)[0]);
}

template <class P>
void ExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                     const TRGSWNTT<P> &trgswntt)
{
    PolynomialNTT<P> decpolyntt;
    TRLWENTT<P> restrlwentt;
    {
        DecomposedNoncePolynomial<P> decpoly;
        NonceDecomposition<P>(decpoly, trlwe[0]);
        TwistINTT<P>(decpolyntt, decpoly[0]);
        for (int m = 0; m < P::k + 1; m++)
#ifdef USE_HEXL
            intel::hexl::EltwiseMultMod(
                &(restrlwentt[m][0].value), &(decpolyntt[0].value),
                &(trgswntt[0][m][0].value), P::n, lvl1P, 1);
#else
            for (int i = 0; i < P::n; i++)
                restrlwentt[m][i] = decpolyntt[i] * trgswntt[0][m][i];
#endif
        for (int i = 1; i < P::l; i++) {
            TwistINTT<P>(decpolyntt, decpoly[i]);
            for (int m = 0; m < P::k + 1; m++)
#ifdef USE_HEXL
            {
                std::array<uint64_t, TFHEpp::lvl1param::n> temp;
                intel::hexl::EltwiseMultMod(temp.data(), &(decpolyntt[0].value),
                                            &(trgswntt[i][m][0].value), P::n,
                                            lvl1P, 1);
                intel::hexl::EltwiseAddMod(&(restrlwentt[m][0].value),
                                           &(restrlwentt[m][0].value),
                                           temp.data(), P::n, lvl1P);
            }
#else
                for (int j = 0; j < P::n; j++)
                    restrlwentt[m][j] += decpolyntt[j] * trgswntt[i][m][j];
#endif
        }
        for (int k = 1; k < P::k; k++) {
            NonceDecomposition<P>(decpoly, trlwe[k]);
            for (int i = 0; i < P::lₐ; i++) {
                TwistINTT<P>(decpolyntt, decpoly[i]);
                for (int m = 0; m < P::k + 1; m++)
#ifdef USE_HEXL
                {
                    std::array<uint64_t, TFHEpp::lvl1param::n> temp;
                    intel::hexl::EltwiseMultMod(
                        temp.data(), &(decpolyntt[0].value),
                        &(trgswntt[i + k * P::l][m][0].value), P::n, lvl1P, 1);
                    intel::hexl::EltwiseAddMod(&(restrlwentt[m][0].value),
                                               &(restrlwentt[m][0].value),
                                               temp.data(), P::n, lvl1P);
                }
#else
                    for (int j = 0; j < P::n; j++)
                        restrlwentt[m][j] +=
                            decpolyntt[j] * trgswntt[i + k * P::lₐ][m][j];
#endif
            }
        }
    }
    DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, trlwe[P::k]);
    for (int i = 0; i < P::l; i++) {
        TwistINTT<P>(decpolyntt, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
#ifdef USE_HEXL
        {
            std::array<uint64_t, TFHEpp::lvl1param::n> temp;
            intel::hexl::EltwiseMultMod(temp.data(), &(decpolyntt[0].value),
                                        &(trgswntt[i + k * P::l][m][0].value),
                                        P::n, lvl1P, 1);
            intel::hexl::EltwiseAddMod(&(restrlwentt[m][0].value),
                                       &(restrlwentt[m][0].value), temp.data(),
                                       P::n, lvl1P);
        }
#else
            for (int j = 0; j < P::n; j++)
                restrlwentt[m][j] +=
                    decpolyntt[j] * trgswntt[i + P::k * P::lₐ][m][j];
#endif
    }
    for (int k = 0; k < P::k + 1; k++) TwistNTT<P>(res[k], restrlwentt[k]);
}

template <class P>
TRGSWFFT<P> ApplyFFT2trgsw(const TRGSW<P> &trgsw)
{
    alignas(64) TRGSWFFT<P> trgswfft;
    for (int i = 0; i < P::k * P::lₐ * P::l̅ₐ + P::l * P::l̅; i++)
        for (int j = 0; j < (P::k + 1); j++)
            TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}

template <class P>
void ApplyFFT2trgsw(TRGSWFFT<P> &trgswfft, const TRGSW<P> &trgsw)
{
    for (int i = 0; i < P::k * P::lₐ * P::l̅ₐ + P::l * P::l̅; i++)
        for (int j = 0; j < (P::k + 1); j++)
            TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
}

template <class P>
HalfTRGSWFFT<P> ApplyFFT2halftrgsw(const HalfTRGSW<P> &trgsw)
{
    alignas(64) HalfTRGSWFFT<P> halftrgswfft;
    for (int i = 0; i < P::l * P::l̅; i++)
        for (int j = 0; j < (P::k + 1); j++)
            TwistIFFT<P>(halftrgswfft[i][j], trgsw[i][j]);
    return halftrgswfft;
}

template <class P>
TRGSWNTT<P> ApplyNTT2trgsw(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < P::k * P::lₐ * P::l̅ₐ + P::l * P::l̅; i++)
        for (int j = 0; j < P::k + 1; j++)
            TwistINTT<P>(trgswntt[i][j], trgsw[i][j]);
    return trgswntt;
}

template <class P>
TRGSWRAINTT<P> ApplyRAINTT2trgsw(const TRGSW<P> &trgsw)
{
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    TRGSWRAINTT<P> trgswntt;
    for (int i = 0; i < P::k * P::lₐ * P::l̅ₐ + P::l * P::l̅; i++)
        for (int j = 0; j < P::k + 1; j++) {
            raintt::TwistINTT<typename P::T, P::nbit, true>(
                trgswntt[i][j], trgsw[i][j], (*raintttable)[1],
                (*raintttwist)[1]);
            for (int k = 0; k < P::n; k++)
                if ((k & ((1 << remainder) - 1)) > 1)
                    trgswntt[i][j][k] =
                        raintt::MulSREDC(trgswntt[i][j][k], raintt::R4);
                else
                    trgswntt[i][j][k] =
                        raintt::MulSREDC(trgswntt[i][j][k], raintt::R2);
        }
    return trgswntt;
}

template <class P>
TRGSWNTT<P> TRGSW2NTT(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < P::k * P::lₐ * P::l̅ₐ + P::l * P::l̅; i++)
        for (int j = 0; j < P::k + 1; j++) {
            PolynomialNTT<P> temp;
            TwistINTT<P>(temp, trgsw[i][j]);
            for (uint32_t k = 0; k < P::n; k++)
                trgswntt[i][j][k] = temp[cuHEpp::BitReverse<P::nbit>(k)];
        }
    return trgswntt;
}

template <class P>
constexpr std::array<typename P::T, P::l> hgen()
{
    std::array<typename P::T, P::l> h{};
    if constexpr (hasq<P>)
        for (int i = 0; i < P::lₐ; i++)
            h[i] = (P::q + (1ULL << ((i + 1) * P::Bgbit - 1))) >>
                   ((i + 1) * P::Bgbit);
    else
        for (int i = 0; i < P::l; i++)
            h[i] = static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                            (i + 1) * P::Bgbit);
    return h;
}

template <class P>
constexpr std::array<typename P::T, P::lₐ> noncehgen()
{
    std::array<typename P::T, P::lₐ> h{};
    if constexpr (hasq<P>)
        for (int i = 0; i < P::lₐ; i++)
            h[i] = (P::q + (1ULL << ((i + 1) * P::Bgₐbit - 1))) >>
                   ((i + 1) * P::Bgₐbit);
    else
        for (int i = 0; i < P::lₐ; i++)
            h[i] = static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                            (i + 1) * P::Bgₐbit);
    return h;
}

// Auxiliary h generation for Double Decomposition (bivariate representation)
// h̅[j] values are used to construct gadget values h[i] * h̅[j] = 2^(width - (i+1)*Bgbit - j*B̅gbit)
// For j=0: no auxiliary shift, so h̅[0] = 1
// For j>0: h̅[j] = 2^(width - j*B̅gbit) which combines with h[i] via modular multiplication
template <class P>
constexpr std::array<typename P::T, P::l̅> h̅gen()
{
    std::array<typename P::T, P::l̅> h̅{};
    h̅[0] = 1;  // j=0 means no auxiliary shift
    for (int i = 1; i < P::l̅; i++)
        h̅[i] = static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                        i * P::B̅gbit);
    return h̅;
}

// Auxiliary h generation for nonce part of TRGSW with Double Decomposition
template <class P>
constexpr std::array<typename P::T, P::l̅ₐ> nonceh̅gen()
{
    std::array<typename P::T, P::l̅ₐ> h̅{};
    h̅[0] = 1;  // j=0 means no auxiliary shift
    for (int i = 1; i < P::l̅ₐ; i++)
        h̅[i] = static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits -
                        i * P::B̅gₐbit);
    return h̅;
}

template <class P>
inline void halftrgswhadd(HalfTRGSW<P> &halftrgsw, const Polynomial<P> &p)
{
    constexpr std::array<typename P::T, P::l> h = hgen<P>();
    constexpr std::array<typename P::T, P::l̅> h̅ = h̅gen<P>();
    for (int i = 0; i < P::l; i++) {
        for (int ī = 0; ī < P::l̅; ī++) {
            for (int j = 0; j < P::n; j++) {
                halftrgsw[i * P::l̅ + ī][P::k][j] +=
                    static_cast<typename P::T>(p[j]) * h[i] * h̅[ī];
            }
        }
    }
}

template <class P>
inline void trgswhadd(TRGSW<P> &trgsw, const Polynomial<P> &p)
{
    constexpr std::array<typename P::T, P::lₐ> nonceh = noncehgen<P>();
    constexpr std::array<typename P::T, P::l̅ₐ> nonceh̅ = nonceh̅gen<P>();
    for (int i = 0; i < P::lₐ; i++) {
        for (int ī = 0; ī < P::l̅ₐ; ī++) {
            for (int k = 0; k < P::k; k++) {
                for (int j = 0; j < P::n; j++) {
                    trgsw[(i * P::l̅ₐ + ī) + k * P::lₐ * P::l̅ₐ][k][j] +=
                        static_cast<typename P::T>(p[j]) * nonceh[i] *
                        nonceh̅[ī];
                }
            }
        }
    }
    constexpr std::array<typename P::T, P::l> h = hgen<P>();
    constexpr std::array<typename P::T, P::l̅> h̅ = h̅gen<P>();
    for (int i = 0; i < P::l; i++) {
        for (int ī = 0; ī < P::l̅; ī++) {
            for (int j = 0; j < P::n; j++) {
                trgsw[(i * P::l̅ + ī) + P::k * P::lₐ * P::l̅ₐ][P::k][j] +=
                    static_cast<typename P::T>(p[j]) * h[i] * h̅[ī];
            }
        }
    }
}

template <class P>
inline void trgswhoneadd(TRGSW<P> &trgsw)
{
    constexpr std::array<typename P::T, P::lₐ> nonceh = noncehgen<P>();
    constexpr std::array<typename P::T, P::l̅ₐ> nonceh̅ = nonceh̅gen<P>();
    for (int i = 0; i < P::lₐ; i++)
        for (int ī = 0; ī < P::l̅ₐ; ī++)
            for (int k = 0; k < P::k; k++)
                trgsw[(i * P::l̅ₐ + ī) + k * P::lₐ * P::l̅ₐ][k][0] +=
                    nonceh[i] * nonceh̅[ī];

    constexpr std::array<typename P::T, P::l> h = hgen<P>();
    constexpr std::array<typename P::T, P::l̅> h̅ = h̅gen<P>();
    for (int i = 0; i < P::l; i++)
        for (int ī = 0; ī < P::l̅; ī++)
            trgsw[(i * P::l̅ + ī) + P::k * P::lₐ * P::l̅ₐ][P::k][0] +=
                h[i] * h̅[ī];
}

template <class P>
void trgswSymEncrypt(TRGSW<P> &trgsw, const Polynomial<P> &p, const double α,
                     const Key<P> &key)
{
    for (TRLWE<P> &trlwe : trgsw) trlweSymEncryptZero<P>(trlwe, α, key);
    trgswhadd<P>(trgsw, p);
}

template <class P>
void trgswSymEncrypt(TRGSW<P> &trgsw, const Polynomial<P> &p, const uint η,
                     const Key<P> &key)
{
    for (TRLWE<P> &trlwe : trgsw) trlweSymEncryptZero<P>(trlwe, η, key);
    trgswhadd<P>(trgsw, p);
}

template <class P>
void trgswSymEncrypt(TRGSW<P> &trgsw, const Polynomial<P> &p,
                     const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trgswSymEncrypt<P>(trgsw, p, P::α, key);
    else
        trgswSymEncrypt<P>(trgsw, p, P::η, key);
}

template <class P>
void halftrgswSymEncrypt(HalfTRGSW<P> &halftrgsw, const Polynomial<P> &p,
                         const double α, const Key<P> &key)
{
    for (TRLWE<P> &trlwe : halftrgsw) trlweSymEncryptZero<P>(trlwe, α, key);
    halftrgswhadd<P>(halftrgsw, p);
}

template <class P>
void halftrgswSymEncrypt(HalfTRGSW<P> &halftrgsw, const Polynomial<P> &p,
                         const uint η, const Key<P> &key)
{
    for (TRLWE<P> &trlwe : halftrgsw) trlweSymEncryptZero<P>(trlwe, η, key);
    halftrgswhadd<P>(halftrgsw, p);
}

template <class P>
void halftrgswSymEncrypt(HalfTRGSW<P> &halftrgsw, const Polynomial<P> &p,
                         const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        halftrgswSymEncrypt<P>(halftrgsw, p, P::α, key);
    else
        halftrgswSymEncrypt<P>(halftrgsw, p, P::η, key);
}

template <class P>
void trgswSymEncrypt(TRGSWFFT<P> &trgswfft, const Polynomial<P> &p,
                     const double α, const Key<P> &key)
{
    TRGSW<P> trgsw;
    trgswSymEncrypt<P>(trgsw, p, α, key);
    ApplyFFT2trgsw<P>(trgswfft, trgsw);
}

template <class P>
void trgswSymEncrypt(TRGSWFFT<P> &trgswfft, const Polynomial<P> &p,
                     const uint η, const Key<P> &key)
{
    TRGSW<P> trgsw;
    trgswSymEncrypt<P>(trgsw, p, η, key);
    ApplyFFT2trgsw<P>(trgswfft, trgsw);
}

template <class P>
void trgswSymEncrypt(TRGSWFFT<P> &trgswfft, const Polynomial<P> &p,
                     const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trgswSymEncrypt<P>(trgswfft, p, P::α, key);
    else
        trgswSymEncrypt<P>(trgswfft, p, P::η, key);
}

template <class P>
void halftrgswSymEncrypt(HalfTRGSWFFT<P> &halftrgswfft, const Polynomial<P> &p,
                         const double α, const Key<P> &key)
{
    HalfTRGSW<P> halftrgsw;
    halftrgswSymEncrypt<P>(halftrgsw, p, α, key);
    halftrgswfft = ApplyFFT2halftrgsw<P>(halftrgsw);
}

template <class P>
void halftrgswSymEncrypt(HalfTRGSWFFT<P> &halftrgswfft, const Polynomial<P> &p,
                         const uint η, const Key<P> &key)
{
    HalfTRGSW<P> halftrgsw;
    halftrgswSymEncrypt<P>(halftrgsw, p, η, key);
    halftrgswfft = ApplyFFT2halftrgsw<P>(halftrgsw);
}

template <class P>
void halftrgswSymEncrypt(HalfTRGSWFFT<P> &halftrgswfft, const Polynomial<P> &p,
                         const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        halftrgswSymEncrypt<P>(halftrgswfft, p, P::α, key);
    else
        halftrgswSymEncrypt<P>(halftrgswfft, p, P::η, key);
}

template <class P>
void trgswSymEncrypt(TRGSWNTT<P> &trgswntt, const Polynomial<P> &p,
                     const double α, const Key<P> &key)
{
    TRGSW<P> trgsw;
    trgswSymEncrypt<P>(trgsw, p, α, key);
    trgswntt = ApplyNTT2trgsw<P>(trgsw);
}

template <class P>
void trgswSymEncrypt(TRGSWNTT<P> &trgswntt, const Polynomial<P> &p,
                     const uint η, const Key<P> &key)
{
    TRGSW<P> trgsw;
    trgswSymEncrypt<P>(trgsw, p, η, key);
    trgswntt = ApplyNTT2trgsw<P>(trgsw);
}

template <class P>
void trgswSymEncrypt(TRGSWNTT<P> &trgswntt, const Polynomial<P> &p,
                     const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trgswSymEncrypt<P>(trgswntt, p, P::α, key);
    else
        trgswSymEncrypt<P>(trgswntt, p, P::η, key);
}

template <class P>
void trgswSymEncrypt(TRGSWRAINTT<P> &trgswraintt, const Polynomial<P> &p,
                     const double α, const Key<P> &key)
{
    TRGSW<P> trgsw;
    trgswSymEncrypt<P>(trgsw, p, α, key);
    trgswraintt = ApplyRAINTT2trgsw<P>(trgsw);
}

template <class P>
void trgswSymEncrypt(TRGSWRAINTT<P> &trgswraintt, const Polynomial<P> &p,
                     const uint η, const Key<P> &key)
{
    if constexpr (hasq<P> && P::q == raintt::P) {
        constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
        constexpr std::array<typename P::T, P::lₐ> h = hgen<P>();
        for (TRLWERAINTT<P> &trlweraintt : trgswraintt) {
            trlweSymEncryptZero<P>(trlweraintt, η, key);
            for (int k = 0; k <= P::k; k++)
                for (int j = 0; j < P::n; j++)
                    trlweraintt[k][j] = raintt::MulSREDC(
                        trlweraintt[k][j], ((j & ((1 << remainder) - 1)) > 1)
                                               ? raintt::R3
                                               : raintt::R2);
        }
        for (int i = 0; i < P::lₐ; i++) {
            Polynomial<P> pscaled;
            for (int j = 0; j < P::n; j++) {
                pscaled[j] =
                    static_cast<uint64_t>(static_cast<int32_t>(p[j]) >= 0
                                              ? p[j]
                                              : (P::q + p[j])) *
                    h[i] % P::q;
            }
            for (int k = 0; k < P::k; k++) {
                PolynomialRAINTT<P> praintt;
                raintt::TwistINTT<typename P::T, P::nbit, false>(
                    praintt, pscaled, (*raintttable)[1], (*raintttwist)[1]);
                for (int j = 0; j < P::n; j++)
                    if ((j & ((1 << remainder) - 1)) > 1)
                        praintt[j] = raintt::MulSREDC(praintt[j], raintt::R4);
                    else
                        praintt[j] = raintt::MulSREDC(praintt[j], raintt::R2);
                for (int j = 0; j < P::n; j++)
                    trgswraintt[i + k * P::lₐ][k][j] = raintt::AddMod(
                        trgswraintt[i + k * P::lₐ][k][j], praintt[j]);
            }
        }
        for (int i = 0; i < P::l; i++) {
            Polynomial<P> pscaled;
            for (int j = 0; j < P::n; j++) {
                pscaled[j] =
                    static_cast<uint64_t>(static_cast<int32_t>(p[j]) >= 0
                                              ? p[j]
                                              : (P::q + p[j])) *
                    h[i] % P::q;
            }
            PolynomialRAINTT<P> praintt;
            raintt::TwistINTT<typename P::T, P::nbit, false>(
                praintt, pscaled, (*raintttable)[1], (*raintttwist)[1]);
            for (int j = 0; j < P::n; j++)
                if ((j & ((1 << remainder) - 1)) > 1)
                    praintt[j] = raintt::MulSREDC(praintt[j], raintt::R4);
                else
                    praintt[j] = raintt::MulSREDC(praintt[j], raintt::R2);
            for (int j = 0; j < P::n; j++)
                trgswraintt[i + P::k * P::lₐ][P::k][j] = raintt::AddMod(
                    trgswraintt[i + P::k * P::lₐ][P::k][j], praintt[j]);
        }
    }
    else {
        TRGSW<P> trgsw = trgswSymEncrypt<P>(p, η, key);
        trgswraintt = ApplyRAINTT2trgsw<P>(trgsw);
    }
}

template <class P>
void trgswSymEncrypt(TRGSWRAINTT<P> &trgswraintt, const Polynomial<P> &p,
                     const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        trgswSymEncrypt<P>(trgswraintt, p, P::α, key);
    else
        trgswSymEncrypt<P>(trgswraintt, p, P::η, key);
}

}  // namespace TFHEpp