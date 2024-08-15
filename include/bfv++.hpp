#pragma once
#include <cstdint>

#include "keyswitch.hpp"
#include "mulfft.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"

// #include "hexl/hexl.hpp"

namespace TFHEpp {

template <class P>
void TRLWEMultWithoutRelinerization(TRLWE3<P> &res, const TRLWE<P> &a,
                                    const TRLWE<P> &b)
{
    PolynomialInFD<P> ffta, fftb, fftc;
    TwistIFFTUInt<P>(ffta, a[0]);
    TwistIFFTUInt<P>(fftb, b[1]);
    MulInFD<P::n>(fftc, ffta, fftb);
    TwistIFFTUInt<P>(ffta, a[1]);
    TwistIFFTUInt<P>(fftb, b[0]);
    FMAInFD<P::n>(fftc, ffta, fftb);
    TwistFFTrescale<P>(res[0], fftc);

    PolyMulRescaleUnsigned<P>(res[1], a[1], b[1]);
    PolyMulRescaleUnsigned<P>(res[2], a[0], b[0]);

    // for (int i = 0; i < P::n; i++) {
    //     uint64_t ri = 0;
    //     for (int j = 0; j <= i; j++)
    //         ri += static_cast<uint64_t>(P::plain_modulus) *
    //               static_cast<uint64_t>(a[0][j]) * b[0][i - j];
    //     for (int j = i + 1; j < P::n; j++)
    //         ri -= P::plain_modulus * static_cast<uint64_t>(a[0][j]) *
    //               b[0][P::n + i - j];
    //     res[2][i] = (ri + (1ULL << 31)) >> 32;
    // }
}

template <class P>
inline void relinKeySwitch(TRLWE<P> &res, const Polynomial<P> &poly,
                           const relinKeyFFT<P> &relinkeyfft)
{
    DecomposedPolynomial<P> decvec;
    Decomposition<P>(decvec, poly);
    PolynomialInFD<P> decvecfft;
    TwistIFFT<P>(decvecfft, decvec[0]);
    TRLWEInFD<P> resfft;
    MulInFD<P::n>(resfft[0], decvecfft, relinkeyfft[0][0]);
    MulInFD<P::n>(resfft[1], decvecfft, relinkeyfft[0][1]);
    for (int i = 1; i < P::l; i++) {
        TwistIFFT<P>(decvecfft, decvec[i]);
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