#pragma once

#include <cstdint>

#include "cloudkey.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"

namespace TFHEpp {

template <class midP, class targetP>
constexpr Polynomial<midP> CBtestvector()
{
    Polynomial<midP> poly;
    constexpr uint32_t bitwidth = bits_needed<targetP::l - 1>();
    for (int i = 0; i < (midP::n >> bitwidth); i++)
        for (int j = 0; j < (1 << bitwidth); j++)
            poly[(i << bitwidth) + j] =
                1ULL << (std::numeric_limits<typename midP::T>::digits -
                         (j + 1) * targetP::Bgbit - 1);
    return poly;
}

template <class P>
constexpr Polynomial<typename P::domainP> CBtestvector()
{
    return CBtestvector<typename P::domainP, typename P::targetP>();
}

template <class bkP, class privksP>
void CircuitBootstrapping(TRGSW<typename privksP::targetP> &trgsw,
                          const TLWE<typename bkP::domainP> &tlwe,
                          const EvalKey &ek)
{
    alignas(64) std::array<TLWE<typename bkP::targetP>, privksP::targetP::l>
        temp;
    GateBootstrappingManyLUT<bkP, privksP::targetP::l>(
        temp, tlwe, ek.getbkfft<bkP>(), CBtestvector<privksP>());
    for (int i = 0; i < privksP::targetP::l; i++) {
        temp[i][privksP::domainP::k * privksP::domainP::n] +=
            1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                     (i + 1) * privksP::targetP::Bgbit - 1);
        for (int k = 0; k < privksP::targetP::k + 1; k++)
            PrivKeySwitch<privksP>(
                trgsw[i + k * privksP::targetP::l], temp[i],
                ek.getprivksk<privksP>("privksk4cb_" + std::to_string(k)));
    }
}

// https://eprint.iacr.org/2024/1318
template <class brP, class ahP>
void AnnihilateCircuitBootstrapping(TRGSW<typename brP::targetP> &trgsw,
                                    const TLWE<typename brP::domainP> &tlwe,
                                    const EvalKey &ek)
{
    static_assert(brP::targetP::k == ahP::k,
                  "brP::targetP::k must be equal to ahP::k");
    alignas(64) std::array<TLWE<typename brP::targetP>, brP::targetP::l> temp;
    GateBootstrappingManyLUT<brP, brP::targetP::l>(
        temp, tlwe, ek.getbkfft<brP>(),
        CBtestvector<typename brP::targetP, typename brP::targetP>());
    for (int i = 0; i < brP::targetP::l; i++) {
        temp[i][brP::targetP::k * brP::targetP::n] +=
            1ULL << (numeric_limits<typename brP::targetP::T>::digits -
                     (i + 1) * brP::targetP::Bgbit - 1);
        TRLWE<typename brP::targetP> temptrlwe;
        InvSampleExtractIndex<typename brP::targetP>(temptrlwe, temp[i], 0);
        AnnihilateKeySwitching<ahP>(
            trgsw[i + brP::targetP::k * brP::targetP::l], temptrlwe,
            ek.getahk<ahP>());
        // Scheme Switching
        for (int k = 0; k < brP::targetP::k; k++)
            trgswfftExternalProduct<ahP>(
                trgsw[i + k * brP::targetP::l],
                trgsw[i + brP::targetP::k * brP::targetP::l],
                ek.getcbsk<ahP>()[k]);
    }
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrapping(TRGSW<typename privksP::targetP> &trgsw,
                          const TLWE<typename iksP::domainP> &tlwe,
                          const EvalKey &ek)
{
    TLWE<typename bkP::domainP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    CircuitBootstrapping<bkP, privksP>(trgsw, tlwelvl0, ek);
}

template <class iksP, class brP, class ahP>
void AnnihilateCircuitBootstrapping(TRGSW<typename brP::targetP> &trgsw,
                                    const TLWE<typename iksP::domainP> &tlwe,
                                    const EvalKey &ek)
{
    TLWE<typename brP::domainP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    AnnihilateCircuitBootstrapping<brP, ahP>(trgsw, tlwelvl0, ek);
}

template <class brP, class privksP>
void CircuitBootstrappingFFT(TRGSWFFT<typename privksP::targetP> &trgswfft,
                             const TLWE<typename brP::domainP> &tlwe,
                             const EvalKey &ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<brP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFT(TRGSWFFT<typename privksP::targetP> &trgswfft,
                             const TLWE<typename iksP::domainP> &tlwe,
                             const EvalKey &ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<iksP, bkP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
}

template <class brP, class ahP>
void AnnihilateCircuitBootstrappingFFT(
    TRGSWFFT<typename brP::targetP> &trgswfft,
    const TLWE<typename brP::domainP> &tlwe, const EvalKey &ek)
{
    alignas(64) TRGSW<typename brP::targetP> trgsw;
    AnnihilateCircuitBootstrapping<brP, ahP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename brP::targetP>(trgswfft, trgsw);
}

template <class iksP, class brP, class ahP>
void AnnihilateCircuitBootstrappingFFT(
    TRGSWFFT<typename brP::targetP> &trgswfft,
    const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
{
    alignas(64) TRGSW<typename brP::targetP> trgsw;
    AnnihilateCircuitBootstrapping<iksP, brP, ahP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename brP::targetP>(trgswfft, trgsw);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingSub(TRGSW<typename privksP::targetP> &trgsw,
                             const TLWE<typename iksP::domainP> &tlwe,
                             const EvalKey &ek)
{
    alignas(64) TLWE<typename bkP::domainP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    alignas(64) std::array<TLWE<typename bkP::targetP>, privksP::targetP::l>
        temp;
    GateBootstrappingManyLUT<bkP, privksP::targetP::l>(
        temp, tlwelvl0, ek.getbkfft<bkP>(), CBtestvector<privksP>());
    for (int i = 0; i < privksP::targetP::l; i++) {
        temp[i][privksP::domainP::k * privksP::domainP::n] +=
            1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                     (i + 1) * privksP::targetP::Bgbit - 1);
        for (int k = 0; k < privksP::targetP::k + 1; k++) {
            alignas(64) TLWE<typename privksP::targetP> subsettlwe;
            SubsetIdentityKeySwitch<privksP>(subsettlwe, temp[i],
                                             ek.getsubiksk<privksP>());
            SubsetPrivKeySwitch<privksP>(
                trgsw[i + k * privksP::targetP::l], subsettlwe,
                ek.getsubprivksk<privksP>("subprivksk4cb_" +
                                          std::to_string(k)));
        }
    }
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingSubFFT(TRGSWFFT<typename privksP::targetP> &trgswfft,
                                const TLWE<typename iksP::domainP> &tlwe,
                                const EvalKey &ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrappingSub<iksP, bkP, privksP>(trgsw, tlwe, ek);
    for (int i = 0; i < (privksP::targetP::k + 1) * privksP::targetP::l; i++)
        for (int j = 0; j < privksP::targetP::k + 1; j++)
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
}

template <class brP, class privksP>
void CircuitBootstrappingFFTInv(
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename brP::domainP> &tlwe, const EvalKey &ek)
{
    alignas(64) TLWE<typename brP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= brP::domainP::k * brP::domainP::n; i++)
        invtlwe[i] = -tlwe[i];
    CircuitBootstrappingFFT<brP, privksP>(invtrgswfft, invtlwe, ek);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFTInv(
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
{
    alignas(64) TLWE<typename iksP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= iksP::domainP::k * iksP::domainP::n; i++)
        invtlwe[i] = -tlwe[i];
    CircuitBootstrappingFFT<iksP, bkP, privksP>(invtrgswfft, invtlwe, ek);
}

template <class brP, class privksP>
void CircuitBootstrappingFFTwithInv(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename brP::domainP> &tlwe, const EvalKey &ek)
{
    constexpr array<typename privksP::targetP::T, privksP::targetP::l> h =
        hgen<typename privksP::targetP>();

    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<brP, privksP>(trgsw, tlwe, ek);
    for (int i = 0; i < (privksP::targetP::k + 1) * privksP::targetP::l; i++)
        for (int j = 0; j < privksP::targetP::k + 1; j++) {
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
            for (int k = 0; k < privksP::targetP::n; k++) trgsw[i][j][k] *= -1;
        }
    for (int i = 0; i < privksP::targetP::l; i++) {
        trgsw[i][0][0] += h[i];
        trgsw[i + privksP::targetP::l][1][0] += h[i];
    }
    for (int i = 0; i < (privksP::targetP::k + 1) * privksP::targetP::l; i++)
        for (int j = 0; j < privksP::targetP::k + 1; j++)
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i][j],
                                                 trgsw[i][j]);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFTwithInv(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<iksP, bkP, privksP>(trgsw, tlwe, ek);
    for (int i = 0; i < (privksP::targetP::k + 1) * privksP::targetP::l; i++)
        for (int j = 0; j < privksP::targetP::k + 1; j++) {
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
            for (int k = 0; k < privksP::targetP::n; k++) trgsw[i][j][k] *= -1;
        }
    trgswhoneadd<typename privksP::targetP>(trgsw);
    for (int i = 0; i < (privksP::targetP::k + 1) * privksP::targetP::l; i++)
        for (int j = 0; j < privksP::targetP::k + 1; j++)
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i][j],
                                                 trgsw[i][j]);
}

}  // namespace TFHEpp