#include "circuitbootstrapping.hpp"

#include <bits/stdint-uintn.h>

namespace TFHEpp {

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingPartial(TRLWE<typename privksP::targetP> &trgswupper,
                                 TRLWE<typename privksP::targetP> &trgswlower,
                                 const TLWE<typename iksP::domainP> &tlwe,
                                 const CircuitKey<bkP, privksP> &ck,
                                 const KeySwitchingKey<iksP> &iksk,
                                 const uint32_t digit)
{
    TLWE<typename bkP::domainP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, iksk);
    const typename bkP::targetP::T μs2 =
        1ULL << (std::numeric_limits<typename privksP::domainP::T>::digits -
                 (digit + 1) * privksP::targetP::Bgbit - 1);
    Polynomial<typename bkP::targetP> testvec;
    for (int i = 0; i < bkP::targetP::n; i++) testvec[i] = μs2;
    TLWE<typename bkP::targetP> tlwemiddle;
    GateBootstrappingTLWE2TLWEFFT<bkP>(tlwemiddle, tlwelvl0, ck.bkfft, testvec);
    tlwemiddle[bkP::targetP::n] += μs2;
    PrivKeySwitch<privksP>(trgswupper, tlwemiddle, ck.privksk[0]);
    PrivKeySwitch<privksP>(trgswlower, tlwemiddle, ck.privksk[1]);
}
#define INST(iksP, bkP, privksP)                                   \
    template void CircuitBootstrappingPartial<iksP, bkP, privksP>( \
        TRLWE<typename privksP::targetP> & trgswupper,       \
        TRLWE<typename privksP::targetP> & trgswlower,       \
        const TLWE<typename iksP::domainP> &tlwe,            \
        const CircuitKey<bkP, privksP> &ck,                  \
        const KeySwitchingKey<iksP> &iksk,                   \
        const uint32_t digit)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

template <class P>
constexpr Polynomial<typename P::domainP> CBtestvector()
{
    Polynomial<typename P::domainP> poly;
    constexpr uint32_t bitwidth = bits_needed<P::targetP::l - 1>();
    for (int i = 0; i < (P::domainP::n >> bitwidth); i++)
        for (int j = 0; j < (1 << bitwidth); j++)
            poly[(i << bitwidth) + j] =
                1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                         (j + 1) * P::targetP::Bgbit - 1);
    return poly;
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrapping(TRGSW<typename privksP::targetP> &trgsw,
                          const TLWE<typename iksP::domainP> &tlwe,
                          const CircuitKey<bkP, privksP> &ck,
                          const KeySwitchingKey<iksP> &iksk)
{
    TLWE<typename bkP::domainP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, iksk);
    std::array<TLWE<typename bkP::targetP>, privksP::targetP::l> temp;
    GateBootstrappingManyLUT<bkP, privksP::targetP::l>(temp, tlwelvl0, ck.bkfft,
                                                       CBtestvector<privksP>());
    for (int i = 0; i < privksP::targetP::l; i++) {
        temp[i][privksP::domainP::n] +=
            1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                     (i + 1) * privksP::targetP::Bgbit - 1);
        PrivKeySwitch<privksP>(trgsw[i], temp[i], ck.privksk[0]);
        PrivKeySwitch<privksP>(trgsw[i + privksP::targetP::l], temp[i],
                               ck.privksk[1]);
    }
}
#define INST(iksP, bkP, privksP)                      \
    template void CircuitBootstrapping<iksP, bkP, privksP>( \
        TRGSW<typename privksP::targetP> & trgsw,     \
        const TLWE<typename iksP::domainP> &tlwe,      \
        const CircuitKey<bkP, privksP> &ck,           \
        const KeySwitchingKey<iksP> &iksk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFT(TRGSWFFT<typename privksP::targetP> &trgswfft,
                             const TLWE<typename iksP::domainP> &tlwe,
                             const CircuitKey<bkP, privksP> &ck,
                             const KeySwitchingKey<iksP> &iksk)
{
    TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<iksP, bkP, privksP>(trgsw, tlwe, ck,iksk);
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
}
#define INST(iksP, bkP, privksP)                               \
    template void CircuitBootstrappingFFT<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe,              \
        const CircuitKey<bkP, privksP> &ck,                    \
        const KeySwitchingKey<iksP> &iksk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFTInv(
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename iksP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck,const KeySwitchingKey<iksP> &iksk)
{
    TLWE<typename iksP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= bkP::domainP::n; i++) invtlwe[i] = -tlwe[i];
    CircuitBootstrappingFFT<iksP,bkP,privksP>(invtrgswfft, invtlwe, ck,iksk);
}
#define INST(iksP, bkP, privksP)                                  \
    template void CircuitBootstrappingFFTInv<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,  \
        const TLWE<typename iksP::domainP> &tlwe,            \
        const CircuitKey<bkP, privksP> &ck,                 \
        const KeySwitchingKey<iksP> &iksk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFTwithInvPartial(
    TRLWEInFD<typename privksP::targetP> &trgswfftupper,
    TRLWEInFD<typename privksP::targetP> &trgswfftlower,
    TRLWEInFD<typename privksP::targetP> &invtrgswfftupper,
    TRLWEInFD<typename privksP::targetP> &invtrgswfftlower,
    const TLWE<typename iksP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck,
    const KeySwitchingKey<iksP> &iksk,
    const uint32_t digit)
{
    constexpr std::array<typename privksP::targetP::T, privksP::targetP::l> h =
        hgen<typename privksP::targetP>();
    TRLWE<typename privksP::targetP> trgswupper, trgswlower;
    CircuitBootstrappingPartial<iksP, bkP, privksP>(trgswupper, trgswlower, tlwe, ck, iksk, digit);
    for (int j = 0; j < 2; j++) {
        TwistIFFT<typename privksP::targetP>(trgswfftupper[j], trgswupper[j]);
        TwistIFFT<typename privksP::targetP>(trgswfftlower[j], trgswlower[j]);
    }
    for (int j = 0; j < privksP::targetP::n; j++) {
        trgswupper[0][j] *= -1;
        trgswupper[1][j] *= -1;
        trgswlower[0][j] *= -1;
        trgswlower[1][j] *= -1;
    }
    trgswupper[0][0] += h[digit];
    trgswlower[1][0] += h[digit];
    for (int j = 0; j < 2; j++) {
        TwistIFFT<typename privksP::targetP>(invtrgswfftupper[j],
                                             trgswupper[j]);
        TwistIFFT<typename privksP::targetP>(invtrgswfftlower[j],
                                             trgswlower[j]);
    }
}
#define INST(iksP, bkP, privksP)                                       \
    template void CircuitBootstrappingFFTwithInvPartial<iksP, bkP, privksP>( \
        TRLWEInFD<typename privksP::targetP> & trgswfftupper,          \
        TRLWEInFD<typename privksP::targetP> & trgswfftlower,          \
        TRLWEInFD<typename privksP::targetP> & invtrgswfftupper,       \
        TRLWEInFD<typename privksP::targetP> & invtrgswfftlower,       \
        const TLWE<typename iksP::domainP> &tlwe,                      \
        const CircuitKey<bkP, privksP> &ck,                            \
        const KeySwitchingKey<iksP> &iksk,                             \
         const uint32_t digit)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingFFTwithInv(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename iksP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck,
    const KeySwitchingKey<iksP> &iksk)
{
    constexpr array<typename privksP::targetP::T, privksP::targetP::l> h =
        hgen<typename privksP::targetP>();

    TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<iksP, bkP, privksP>(trgsw, tlwe, ck, iksk);
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++) {
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
            for (int k = 0; k < privksP::targetP::n; k++) trgsw[i][j][k] *= -1;
        }
    for (int i = 0; i < privksP::targetP::l; i++) {
        trgsw[i][0][0] += h[i];
        trgsw[i + privksP::targetP::l][1][0] += h[i];
    }
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i][j],
                                                 trgsw[i][j]);
}
#define INST(iksP, bkP, privksP)                                      \
    template void CircuitBootstrappingFFTwithInv<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,         \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,      \
        const TLWE<typename iksP::domainP> &tlwe,               \
        const CircuitKey<bkP, privksP> &ck,                     \
        const KeySwitchingKey<iksP> &iksk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

}  // namespace TFHEpp
