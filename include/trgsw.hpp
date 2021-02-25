#pragma once

#include <array>
#include <cstdint>
#include <params.hpp>
#include <mulfft.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
using namespace std;

template<class P>
constexpr typename P::T offsetgen(){
    typename P::T offset = 0;
    for (int i = 1; i <= P::l; i++)
        offset += P::Bg / 2 * (1ULL << (numeric_limits<typename P::T>::digits - i * P::Bgbit));
    return offset;
}

template <class P>
inline void Decomposition(DecomposedTRLWE<P> &decvec,
                          const TRLWE<P> &trlwe)
{
    constexpr typename P::T offset = offsetgen<P>();
    constexpr typename P::T mask = static_cast<typename P::T>((1ULL << P::Bgbit) - 1);
    constexpr typename P::T halfBg = (1ULL << (P::Bgbit - 1));

    for (int j = 0; j < P::n; j++) {
        typename P::T temp0 = trlwe[0][j] + offset;
        typename P::T temp1 = trlwe[1][j] + offset;
        for (int i = 0; i < P::l; i++)
            decvec[i][j] =
                ((temp0 >> (numeric_limits<typename P::T>::digits - (i + 1) * P::Bgbit)) &
                 mask) -
                halfBg;
        for (int i = 0; i < P::l; i++)
            decvec[i + P::l][j] =
                ((temp1 >> (numeric_limits<typename P::T>::digits - (i + 1) * P::Bgbit)) &
                 mask) -
                halfBg;
    }
}

template <class P>
inline void DecompositionFFT(DecomposedTRLWEInFD<P> &decvecfft,
                                 const TRLWE<P> &trlwe)
{
    DecomposedTRLWE<P> decvec;
    Decomposition<P>(decvec, trlwe);
    for (int i = 0; i < 2 * P::l; i++) TwistIFFT<P>(decvecfft[i], decvec[i]);
}

template<class P>
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

template<class P>
inline TRGSW<P> trgswSymEncrypt(const typename make_signed<typename P::T>::type p, const double α, const Key<P> &key)
{
    array<typename P::T, P::l> h;
    for (int i = 0; i < P::l; i++) h[i] = 1ULL << (numeric_limits<typename P::T>::digits - (i + 1) * P::Bgbit);

    TRGSW<P> trgsw;
    for (TRLWE<P> &trlwe : trgsw) trlwe = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::l; i++) {
        trgsw[i][0][0] += static_cast<typename P::T>(p) * h[i];
        trgsw[i + P::l][1][0] += static_cast<typename P::T>(p) * h[i];
    }
    return trgsw;
}

template<class P>
TRGSWFFT<P> trgswfftSymEncrypt(const typename make_signed<typename P::T>::type p, const double α, const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, α, key);
    TRGSWFFT<P> trgswfft;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}

TRGSWFFT<lvl1param> trgswfftSymEncryptlvl1(const typename make_signed<lvl1param::T>::type p, const double α, const Key<lvl1param> &key);
void trgswfftExternalProductlvl1(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &trlwe,
                                 const TRGSWFFT<lvl1param> &trgswfft);

TRGSWFFT<lvl2param> trgswfftSymEncryptlvl2(const typename make_signed<lvl2param::T>::type p, const double α, const Key<lvl2param> &key);
void trgswfftExternalProductlvl2(TRLWE<lvl2param> &res, const TRLWE<lvl2param> &trlwe,
                                 const TRGSWFFT<lvl2param> &trgswfft);
}  // namespace TFHEpp