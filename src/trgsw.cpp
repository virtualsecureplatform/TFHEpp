#include <array>
#include <limits>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace SPQLIOSpp;

inline TRGSW<lvl1param> trgswSymEncryptlvl1(
    const make_signed<lvl1param::T>::type p, const double α,
    const Key<lvl1param> &key)
{
    return trgswSymEncrypt<lvl1param>(p, α, key);
}

inline TRGSW<lvl2param> trgswSymEncryptlvl2(
    const make_signed<lvl2param::T>::type p, const double α,
    const Key<lvl2param> &key)
{
    return trgswSymEncrypt<lvl2param>(p, α, key);
}

TRGSWFFT<lvl1param> trgswfftSymEncryptlvl1(
    const make_signed<lvl1param::T>::type p, const double α,
    const Key<lvl1param> &key)
{
    return trgswfftSymEncrypt<lvl1param>(p, α, key);
}

TRGSWFFT<lvl2param> trgswfftSymEncryptlvl2(
    const make_signed<lvl2param::T>::type p, const double α,
    const Key<lvl2param> &key)
{
    return trgswfftSymEncrypt<lvl2param>(p, α, key);
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
    DecompositionFFT<lvl1param>(decvecfft, trlwe);
}

inline void DecompositionFFTlvl2(DecomposedTRLWEInFD<lvl2param> &decvecfft,
                                 const TRLWE<lvl2param> &trlwe)
{
    DecompositionFFT<lvl2param>(decvecfft, trlwe);
}

void trgswfftExternalProductlvl1(TRLWE<lvl1param> &res,
                                 const TRLWE<lvl1param> &trlwe,
                                 const TRGSWFFT<lvl1param> &trgswfft)
{
    trgswfftExternalProduct<lvl1param>(res, trlwe, trgswfft);
}

void trgswfftExternalProductlvl2(TRLWE<lvl2param> &res,
                                 const TRLWE<lvl2param> &trlwe,
                                 const TRGSWFFT<lvl2param> &trgswfft)
{
    trgswfftExternalProduct<lvl2param>(res, trlwe, trgswfft);
}
}  // namespace TFHEpp