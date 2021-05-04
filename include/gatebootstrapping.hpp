#pragma once

#include "./cloudkey.hpp"
#include "./detwfa.hpp"
#include "./utils.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
void RotatedTestVector(array<array<typename P::T, P::n>, 2> &testvector,
                       uint32_t bara, const typename P::T μ);

template <class P>
void GateBootstrappingTLWE2TRLWEFFT(TRLWE<typename P::targetP> &acc,
                                    const TLWE<typename P::domainP> &tlwe,
                                    const BootstrappingKeyFFT<P> &bkfft);

template <class P>
void GateBootstrappingTLWE2TLWEFFT(TLWE<typename P::targetP> &res,
                                   const TLWE<typename P::domainP> &tlwe,
                                   const BootstrappingKeyFFT<P> &bkfft);

template <class P>
void GateBootstrappingTLWE2TLWEFFTvariableMu(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft, const typename P::targetP::T μs2);

void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const GateKey &gk);
}  // namespace TFHEpp