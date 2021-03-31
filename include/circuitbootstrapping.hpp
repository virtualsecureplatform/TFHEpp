#pragma once

#include <cloudkey.hpp>
#include <cstdint>
#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {
using namespace std;

template <class bkP, class privksP>
inline void CircuitBootstrappingPartial(TRLWE<typename privksP::targetP> &trgswupper,
                                 TRLWE<typename privksP::targetP> &trgswlower,
                                 const TLWE<typename bkP::domainP> &tlwe,
                                 const CircuitKey<bkP, privksP> &ck,
                                 const uint32_t digit)
{
    TLWE<typename bkP::targetP> tlwemiddle;
    GateBootstrappingTLWE2TLWEFFTvariableMu<bkP>(
        tlwemiddle, tlwe, ck.bkfft,
        1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                    (digit + 1) * privksP::targetP::Bgbit - 1));
    PrivKeySwitch<privksP>(trgswupper, tlwemiddle, ck.privksk[0]);
    PrivKeySwitch<privksP>(trgswlower, tlwemiddle,
                            ck.privksk[1]);
}

template <class bkP, class privksP>
inline void CircuitBootstrapping(TRGSW<typename privksP::targetP> &trgsw,
                                 const TLWE<typename bkP::domainP> &tlwe,
                                 const CircuitKey<bkP, privksP> &ck)
{
    for (int i = 0; i < privksP::targetP::l; i++) {
        CircuitBootstrappingPartial(trgsw[i],trgsw[i+privksP::targetP::l],tlwe,ck,i);
    }
}

template <class bkP, class privksP>
inline void CircuitBootstrappingFFT(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    const TLWE<typename bkP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck)
{
    for (int i = 0; i < privksP::targetP::l; i++){
        TRLWE<typename privksP::targetP> trgswupper, trgswlower;
        CircuitBootstrappingPartial<bkP, privksP>(trgswupper, trgswlower, tlwe, ck,i);
        for (int j = 0; j < 2; j++){
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgswupper[j]);
            TwistIFFT<typename privksP::targetP>(trgswfft[i+privksP::targetP::l][j], trgswlower[j]);
        }
    }
}

template <class bkP, class privksP>
inline void CircuitBootstrappingFFTInv(
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename bkP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck)
{
    TLWE<typename bkP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= bkP::domainP::n; i++) invtlwe[i] = -tlwe[i];
    CircuitBootstrappingFFT(invtrgswfft, tlwe, ck);
}

template <class bkP, class privksP>
inline void CircuitBootstrappingFFTwithInv(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename bkP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck)
{
    constexpr array<typename privksP::targetP::T, privksP::targetP::l> h =
        hgen<typename privksP::targetP>();
    for (int i = 0; i < privksP::targetP::l; i++){
        TRLWE<typename privksP::targetP> trgswupper, trgswlower;
        CircuitBootstrappingPartial(trgswupper,trgswlower,tlwe,ck,i);
        for (int j = 0; j < 2; j++){
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgswupper[j]);
            TwistIFFT<typename privksP::targetP>(trgswfft[i+privksP::targetP::l][j], trgswlower[j]);
        }
        for (int j = 0; j < privksP::targetP::n; j++) {
            trgswupper[0][j] *= -1;
            trgswupper[1][j] *= -1;
            trgswlower[0][j] *= -1;
            trgswlower[1][j] *= -1;
        }
        trgswupper[0][0] += h[i];
        trgswlower[1][0] += h[i];
        for (int j = 0; j < 2; j++){
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i][j], trgswupper[j]);
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i+privksP::targetP::l][j], trgswlower[j]);
        }
    }
}

}  // namespace TFHEpp