#pragma once
#include "keyswitch.hpp"
#include "mulfft.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"

namespace TFHEpp {
template <class P>
inline void RemoveSign(UnsignedPolynomial<P> &res, const Polynomial<P> &a)
{
    for (int i = 0; i < P::n; i++) res[i] = (a[i] + 1) / 2;
}
template <class P>
void TRLWEMultWithoutRelinerization(TRLWE3<P> &res, const TRLWE<P> &a,
                                    const TRLWE<P> &b)
{
    UnsignedTRLWE<P> aa, bb;
    RemoveSign<P>(aa[0], a[0]);
    RemoveSign<P>(aa[1], a[1]);
    RemoveSign<P>(bb[0], b[0]);
    RemoveSign<P>(bb[1], b[1]);

    PolynomialInFD<P> ffta, fftb, fftc;
    TwistIFFT<P>(ffta, aa[0]);
    TwistIFFT<P>(fftb, bb[1]);
    MulInFD<P::n>(fftc, ffta, fftb);
    TwistIFFT<P>(ffta, aa[1]);
    TwistIFFT<P>(fftb, bb[0]);
    FMAInFD<P::n>(fftc, ffta, fftb);
    TwistFFTrescale<P>(res[0], fftc);

    PolyMulRescaleUnsigned<P>(res[1], aa[1], bb[1]);

    PolyMulRescaleUnsigned<P>(res[2], aa[0], bb[0]);
}

template <class P>
inline void relinKeySwitch(TRLWE<P> &res, const Polynomial<P> &poly,
                           const relinKeyFFT<P> &relinkeyfft)
{
    DecomposedPolynomialInFD<P> decvecfft;
    DecompositionPolynomialFFT<P>(decvecfft, poly, 0);
    TRLWEInFD<P> resfft;
    MulInFD<P::n>(resfft[0], decvecfft, relinkeyfft[0][0]);
    MulInFD<P::n>(resfft[1], decvecfft, relinkeyfft[0][1]);
    for (int i = 1; i < P::l; i++) {
        DecompositionPolynomialFFT<P>(decvecfft, poly, i);
        FMAInFD<P::n>(resfft[0], decvecfft, relinkeyfft[i][0]);
        FMAInFD<P::n>(resfft[1], decvecfft, relinkeyfft[i][1]);
    }
    TwistFFT<P>(res[0], resfft[0]);
    TwistFFT<P>(res[1], resfft[1]);
}

template <class P>
inline void Relinearization(TRLWE<P> &res, const TRLWE3<P> &mult,
                            const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE<P> squareterm;
    relinKeySwitch<P>(squareterm, mult[2], relinkeyfft);
    for (int i = 0; i < P::n; i++) res[0][i] = mult[0][i] + squareterm[0][i];
    for (int i = 0; i < P::n; i++) res[1][i] = mult[1][i] + squareterm[1][i];
}

template <class P>
inline void TRLWEMult(TRLWE<P> &res, const TRLWE<P> &a, const TRLWE<P> &b,
                      const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE3<P> resmult;
    TRLWEMultWithoutRelinerization<P>(resmult, a, b);
    Relinearization<P>(res, resmult, relinkeyfft);
}

template <class P>
inline void TLWEMult(TLWE<typename P::targetP> &res,
                     const TLWE<typename P::domainP> &a,
                     const TLWE<typename P::domainP> &b,
                     const relinKeyFFT<typename P::targetP> &relinkeyfft,
                     const PrivateKeySwitchingKey<P> &privksk)
{
    TRLWE<typename P::targetP> trlweres, trlwea, trlweb;
    PrivKeySwitch<P>(trlwea, a, privksk);
    PrivKeySwitch<P>(trlweb, b, privksk);
    TRLWEMult<typename P::targetP>(trlweres, trlwea, trlweb, relinkeyfft);
    SampleExtractIndex<typename P::targetP>(res, trlweres, 0);
}
}  // namespace TFHEpp