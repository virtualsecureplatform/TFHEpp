#pragma once

#include <cloudkey.hpp>
#include <detwfa.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

template <class P>
inline void RotatedTestVector(array<array<typename P::T, P::n>, 2> &testvector,
                              uint32_t bara, const typename P::T μ)
{
    testvector[0] = {};
    if (bara < P::n) {
        for (int i = 0; i < bara; i++) testvector[1][i] = -μ;
        for (int i = bara; i < P::n; i++) testvector[1][i] = μ;
    }
    else {
        const typename P::T baraa = bara - P::n;
        for (int i = 0; i < baraa; i++) testvector[1][i] = μ;
        for (int i = baraa; i < P::n; i++) testvector[1][i] = -μ;
    }
}

void GateBootstrappingTLWE2TRLWEFFTlvl01(
    TRLWE<lvl1param> &acc, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl01param> &bkfft);
void GateBootstrappingTLWE2TLWEFFTlvl01(
    TLWE<lvl1param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl01param> &bkfft);

template <class P>
inline void GateBootstrappingTLWE2TLWEFFTvariableMu(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootStrappingKeyFFT<P> &bkfft, const typename P::targetP::T μs2)
{
    TRLWE<typename P::targetP> acc, temp;
    uint32_t bara = 2 * P::targetP::n - modSwitchFromTorus<typename P::targetP>(
                                            tlwe[P::domainP::n]);
    RotatedTestVector<typename P::targetP>(acc, bara, μs2);
    for (int i = 0; i < P::domainP::n; i++) {
        bara = modSwitchFromTorus<typename P::targetP>(tlwe[i]);
        if (bara == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOne<typename P::targetP>(
            acc, bkfft[i], bara);
    }
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
    res[P::targetP::n] += μs2;
}
void GateBootstrappingTLWE2TLWEFFTvariableMulvl02(
    TLWE<lvl2param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl02param> &bkfftlvl02, const lvl2param::T μs2);
void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const GateKey &gk);
}  // namespace TFHEpp