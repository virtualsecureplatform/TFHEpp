#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;


template <class P>
inline void PolynomialMulByXaiMinusOne(array<typename P::T, P::n> &res,
                                       const array<typename P::T, P::n> &poly, const typename P::T a)
{
    if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n] - poly[i];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n] - poly[i];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa] - poly[i];
    }
}

void PolynomialMulByXaiMinusOnelvl1(Polynomial<lvl1param> &res,
                                    const Polynomial<lvl1param> &poly,
                                    const lvl1param::T a);

template <class P>
inline void RotatedTestVector(array<array<typename P::T, P::n>, 2> &testvector, uint32_t bara,
                              const typename P::T μ)
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

void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWE<lvl1param> &acc, const TLWE<lvl0param> &tlwe,
                                         const BootStrappingKeyFFT<lvl01param> &bkfft);
void GateBootstrappingTLWE2TLWEFFTlvl01(TLWE<lvl1param> &res, const TLWE<lvl0param> &tlwe,
                                        const BootStrappingKeyFFT<lvl01param> &bkfft);

template<class P>
inline void GateBootstrappingTLWE2TLWEFFTvariableMu(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootStrappingKeyFFT<P> &bkfft, const typename P::targetP::T μs2)
{
    TRLWE<typename P::targetP> acc,temp;
    uint32_t bara =
        2 * P::targetP::n - modSwitchFromTorus<typename P::targetP>(tlwe[P::domainP::n]);
    RotatedTestVector<typename P::targetP>(acc, bara, μs2);
    for (int i = 0; i < P::domainP::n; i++) {
        bara = modSwitchFromTorus<typename P::targetP>(tlwe[i]);
        if (bara == 0) continue;
        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[1], acc[1], bara);
        trgswfftExternalProduct<typename P::targetP>(temp, temp, bkfft[i]);
        for (int i = 0; i < P::targetP::n; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
    res[P::targetP::n] += μs2;
}
void GateBootstrappingTLWE2TLWEFFTvariableMulvl02(
    TLWE<lvl2param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl02param> &bkfftlvl02, const lvl2param::T μs2);
void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe, const GateKey &gk);
}  // namespace TFHEpp