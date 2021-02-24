#include <array>
#include <limits>
#include <mulfft.hpp>
#include <params.hpp>
#include <trlwe.hpp>
#include <utils.hpp>
#include <trgsw.hpp>

namespace TFHEpp {

template<class P>
inline TRGSW<P> trgswSymEncrypt(make_signed<typename P::T> p, double α, Key<P> &key)
{
    array<P::T, P::l> h;
    for (int i = 0; i < P::l; i++) h[i] = 1ULL << (numeric_limits<P::T>::digits - (i + 1) * P::Bgbit);

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