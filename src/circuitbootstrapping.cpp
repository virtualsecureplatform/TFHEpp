#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>

namespace TFHEpp {
using namespace std;

template <class bkP, class privksP>
inline void CircuitBootstrapping(TRGSW<typename privksP::targetP> &trgsw,
                                 const TLWE<typename bkP::domainP> &tlwe,
                                 const CircuitKey<bkP, privksP> &ck)
{
    TLWE<typename bkP::targetP> tlwemiddle;
    for (int i = 0; i < privksP::targetP::l; i++) {
        GateBootstrappingTLWE2TLWEFFTvariableMu<bkP>(
            tlwemiddle, tlwe, ck.bkfft,
            1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                     (i + 1) * privksP::targetP::Bgbit - 1));
        PrivKeySwitch<privksP>(trgsw[i], tlwemiddle, 0, ck.privksk);
        PrivKeySwitch<privksP>(trgsw[i + privksP::targetP::l], tlwemiddle, 1,
                               ck.privksk);
    }
}

template <class bkP, class privksP>
inline void CircuitBootstrappingFFT(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    const TLWE<typename bkP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck)
{
    TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<bkP, privksP>(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
}

void CircuitBootstrappingFFTlvl01(TRGSWFFT<lvl1param> &trgswfft,
                                  const TLWE<lvl0param> &tlwe,
                                  const CircuitKey<lvl02param, lvl21param> &ck)
{
    CircuitBootstrappingFFT<lvl02param, lvl21param>(trgswfft, tlwe, ck);
}

void CircuitBootstrappingFFTlvl02(TRGSWFFT<lvl2param> &trgswfft,
                                  const TLWE<lvl0param> &tlwe,
                                  const CircuitKey<lvl02param, lvl22param> &ck)
{
    CircuitBootstrappingFFT<lvl02param, lvl22param>(trgswfft, tlwe, ck);
}

template <class bkP, class privksP>
inline void CircuitBootstrappingFFTInv(
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename bkP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck)
{
    TLWE<typename bkP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= bkP::domainP::n; i++) invtlwe[i] = -tlwe[i];
    TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<bkP, privksP>(trgsw, invtlwe, ck);
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i][j],
                                                 trgsw[i][j]);
}

void CircuitBootstrappingFFTInvlvl01(
    TRGSWFFT<lvl1param> &invtrgswfft, const TLWE<lvl0param> &tlwe,
    const CircuitKey<lvl02param, lvl21param> &ck)
{
    CircuitBootstrappingFFTInv<lvl02param, lvl21param>(invtrgswfft, tlwe, ck);
}

template <class P>
constexpr array<typename P::T, P::l> hgen()
{
    array<typename P::T, P::l> h;
    for (int i = 0; i < P::l; i++)
        h[i] = 1ULL << (numeric_limits<typename P::T>::digits -
                        (i + 1) * P::Bgbit);
    return h;
}

template <class bkP, class privksP>
inline void CircuitBootstrappingFFTwithInv(
    TRGSWFFT<typename privksP::targetP> &trgswfft,
    TRGSWFFT<typename privksP::targetP> &invtrgswfft,
    const TLWE<typename bkP::domainP> &tlwe, const CircuitKey<bkP, privksP> &ck)
{
    TRGSW<typename privksP::targetP> trgsw;
    array<typename privksP::targetP::T, privksP::targetP::l> h =
        hgen<typename privksP::targetP>();
    CircuitBootstrapping<bkP, privksP>(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<typename privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
    for (int i = 0; i < privksP::targetP::l; i++) {
        for (int j = 0; j < privksP::targetP::n; j++) {
            trgsw[i][0][j] *= -1;
            trgsw[i][1][j] *= -1;
            trgsw[i + privksP::targetP::l][0][j] *= -1;
            trgsw[i + privksP::targetP::l][1][j] *= -1;
        }
        trgsw[i][0][0] += h[i];
        trgsw[i + privksP::targetP::l][1][0] += h[i];
    }
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<typename privksP::targetP>(invtrgswfft[i][j],
                                                 trgsw[i][j]);
}

void CircuitBootstrappingFFTwithInvlvl01(
    TRGSWFFT<lvl1param> &trgswfft, TRGSWFFT<lvl1param> &invtrgswfft,
    const TLWE<lvl0param> &tlwe, const CircuitKey<lvl02param, lvl21param> &ck)
{
    CircuitBootstrappingFFTwithInv<lvl02param, lvl21param>(
        trgswfft, invtrgswfft, tlwe, ck);
}
}  // namespace TFHEpp