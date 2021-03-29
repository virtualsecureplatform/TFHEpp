#pragma once

#include <array>
#include <cstdint>
#include <mulfft.hpp>
#include <params.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
using namespace std;

template <class P>
constexpr typename P::T offsetgen()
{
    typename P::T offset = 0;
    for (int i = 1; i <= P::l; i++)
        offset +=
            P::Bg / 2 *
            (1ULL << (numeric_limits<typename P::T>::digits - i * P::Bgbit));
    return offset;
}


template <class P>
inline void DecompositionPolynomial(DecomposedPolynomial<P> &decpoly, const Polynomial<P> &poly)
{
    constexpr typename P::T offset = offsetgen<P>();
    constexpr typename P::T mask =
        static_cast<typename P::T>((1ULL << P::Bgbit) - 1);
    constexpr typename P::T halfBg = (1ULL << (P::Bgbit - 1));

    for (int j = 0; j < P::n; j++) {
        typename P::T temp = poly[j] + offset;
        for (int i = 0; i < P::l; i++)
            decpoly[i][j] = ((temp >> (numeric_limits<typename P::T>::digits -
                                       (i + 1) * P::Bgbit)) &
                            mask) -
                           halfBg;
    }
}

template <class P>
inline void DecompositionFFT(DecomposedTRLWEInFD<P> &decvecfft,
                             const TRLWE<P> &trlwe)
{
    DecomposedPolynomial<P> decpoly;
    DecompositionPolynomial<P>(decpoly, trlwe[0]);
    for (int i = 0; i < P::l; i++) TwistIFFT<P>(decvecfft[i], decpoly[i]);
    DecompositionPolynomial<P>(decpoly, trlwe[1]);
    for (int i = 0; i < P::l; i++) TwistIFFT<P>(decvecfft[i+P::l], decpoly[i]);
}

template <class P>
void trgswfftExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWFFT<P> &trgswfft)
{
    DecomposedTRLWEInFD<P> decvecfft;
    DecompositionFFT<P>(decvecfft, trlwe);
    TRLWEInFD<P> restrlwefft;
    MulInFD<P::n>(restrlwefft[0], decvecfft[0], trgswfft[0][0]);
    MulInFD<P::n>(restrlwefft[1], decvecfft[0], trgswfft[0][1]);
    for (int i = 1; i < 2 * P::l; i++) {
        FMAInFD<P::n>(restrlwefft[0], decvecfft[i], trgswfft[i][0]);
        FMAInFD<P::n>(restrlwefft[1], decvecfft[i], trgswfft[i][1]);
    }
    TwistFFT<P>(res[0], restrlwefft[0]);
    TwistFFT<P>(res[1], restrlwefft[1]);
}

template <class P>
inline constexpr array<typename P::T, P::l> hgen()
{
    array<typename P::T, P::l> h{};
    for (int i = 0; i < P::l; i++)
        h[i] = 1ULL << (numeric_limits<typename P::T>::digits -
                        (i + 1) * P::Bgbit);
    return h;
}

template <class P>
inline TRGSW<P> trgswSymEncrypt(
    const typename make_signed<typename P::T>::type p, const double α,
    const Key<P> &key)
{
    constexpr array<typename P::T, P::l> h = hgen<P>();

    TRGSW<P> trgsw;
    for (TRLWE<P> &trlwe : trgsw) trlwe = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::l; i++) {
        trgsw[i][0][0] += static_cast<typename P::T>(p) * h[i];
        trgsw[i + P::l][1][0] += static_cast<typename P::T>(p) * h[i];
    }
    return trgsw;
}

template <class P>
TRGSWFFT<P> trgswfftSymEncrypt(
    const typename make_signed<typename P::T>::type p, const double α,
    const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, α, key);
    TRGSWFFT<P> trgswfft;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}
}  // namespace TFHEpp