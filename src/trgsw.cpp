#include <array>
#include <limits>
#include <mulfft.hpp>
#include <params.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {

template<class P>
inline TRGSW<P> trgswSymEncrypt(make_signed<typename P::T> p, double α, Key<P> &key)
{
    array<P::T, P::l> h;
    if constexpr (P::T==uint32_t) for (int i = 0; i < P::l; i++) h[i] = 1U << (32 - (i + 1) * P::Bgbit);
    else if constexpr(P::T==uint64_t) for (int i = 0; i < P::l; i++) h[i] = 1ULL << (64 - (i + 1) * P::Bgbit);
    else static_assert(false_v<T>, "Undefined trgswSymEncrypt!");

    TRGSW<P> trgsw;
    for (TRLWElvl1 &trlwe : trgsw) trlwe = trlweSymEncryptZerolvl1(α, key);
    for (int i = 0; i < P::l; i++) {
        trgsw[i][0][0] += static_cast<P::T>(p) * h[i];
        trgsw[i + P::l][1][0] += static_cast<P::T>(p) * h[i];
    }
    return trgsw;
}

inline TRGSW<lvl1param> trgswSymEncryptlvl1(make_signed<lvl1param::T> p, double α, Key<lvl1param> &key)
{
    return trgswSymEncrypt<lvl1param>(p,α,key);
}

inline TRGSW<lvl2param> trgswSymEncryptlvl2(make_signed<lvl2param::T> p, double α, Key<lvl2param> &key)
{
    return trgswSymEncrypt<lvl2param>(p,α,key);
}

template<class P>
TRGSWFFT<P> trgswfftSymEncrypt(make_signed<typename P::T> p, double α, Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, α, key);
    TRGSWFFT<P> trgswfft;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}

TRGSWFFT<lvl1param> trgswfftSymEncryptlvl1(make_signed<lvl1param::T> p, double α, Key<lvl1param> &key)
{
    return trgswfftSymEncrypt<lvl1param>(p,α,key);
}

TRGSWFFT<lvl2param> trgswfftSymEncryptlvl2(make_signed<lvl2param::T> p, double α, Key<lvl2param> &key)
{
    return trgswfftSymEncrypt<lvl2param>(p,α,key);
}

template<class P>
constexpr typename P::T offsetgen(){
    P::T offset = 0;
    for (int i = 1; i <= P::l; i++)
        offset += P::Bg / 2 * (1ULL << (numeric_limits<P::T>::digits - i * P::Bgbit));
    return offset;
}

template <class P>
inline void Decomposition(DecomposedTRLWE<P> &decvec,
                          const TRLWE<P> &trlwe)
{
    constexpr P::T offset = offsetgen<P>();
    constexpr P::T mask = static_cast<P::T>((1ULL << P::Bgbit) - 1);
    constexpr P::T halfBg = (1ULL << (P::Bgbit - 1));

    for (int j = 0; j < P::n; j++) {
        P::T temp0 = trlwe[0][j] + offset;
        P::T temp1 = trlwe[1][j] + offset;
        for (int i = 0; i < P::l; i++)
            decvec[i][j] =
                ((temp0 >> (numeric_limits<P::T>::digits - (i + 1) * P::Bgbit)) &
                 mask) -
                halfBg;
        for (int i = 0; i < P::l; i++)
            decvec[i + P::l][j] =
                ((temp1 >> (numeric_limits<T>::digits - (i + 1) * P::Bgbit)) &
                 mask) -
                halfBg;
    }
}

inline void Decompositionlvl1(DecomposedTRLWE<lvl1param> &decvec,
                              const TRLWE<lvl1param> &trlwe)
{
    Decomposition<lvl1param>(decvec, trlwe);
}

inline void Decompositionlvl2(DecomposedTRLWE<lvl2param> &decvec,
                              const TRLWE<lvl2param> &trlwe)
{
    Decomposition<lvl2param>(decvec, trlwe);
}

template <class P>
inline void DecompositionFFT(DecomposedTRLWEInFD<P> &decvecfft,
                                 const TRLWE<P> &trlwe)
{
    DecomposedTRLWE<P> decvec;
    Decomposition<P>(decvec, trlwe);
    for (int i = 0; i < 2 * P::l; i++) TwistIFFT<P>(decvecfft[i], decvec[i]);
}

inline void DecompositionFFTlvl1(DecomposedTRLWEInFD<lvl1param> &decvecfft,
                                 const TRLWE<lvl1param> &trlwe)
{
    DecompositionFFT<lvl1param>(decvecfft,trlwe);
}

inline void DecompositionFFTlvl2(DecomposedTRLWEInFD<lvl2param> &decvecfft,
                                 const TRLWE<lvl2param> &trlwe)
{
    DecompositionFFT<lvl2param>(decvecfft,trlwe);
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

void trgswfftExternalProductlvl1(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &trlwe,
                                 const TRGSWFFT<lvl1param> &trgswfft)
{
    trgswfftExternalProduct<lvl1param>(res,trlwe,trgswfft);
}

void trgswfftExternalProductlvl2(TRLWE<lvl2param> &res, const TRLWE<lvl2param> &trlwe,
                                 const TRGSWFFT<lvl2param> &trgswfft)
{
    trgswfftExternalProduct<lvl2param>(res,trlwe,trgswfft);
}
}  // namespace TFHEpp