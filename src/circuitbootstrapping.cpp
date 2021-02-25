#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>

namespace TFHEpp {
using namespace std;

template<class bkP,class privksP>
inline void CircuitBootstrapping(TRGSW<typename privksP::targetP> &trgsw, const TLWE<typename bkP::domainP>  &tlwe,
                                 const CircuitKey<bkP,privksP> &ck)
{
    TLWE<typename bkP::targetP> tlwe;
    for (int i = 0; i < privksP::targetP::l; i++) {
        GateBootstrappingTLWE2TLWEFFT<bkP>(
            tlwe, tlwe, ck.bkfft,
            1ULL << (numeric_limits<typename privksP::domainP::T>::digits - (i + 1) * privksP::targetP::Bgbit - 1));
        PrivKeySwitchlvl21(trgsw[i], tlwe, 0, ck.privksk);
        PrivKeySwitchlvl21(trgsw[i + privksP::targetP::l], tlwe, 1, ck.privksk);
    }
}

template<class bkP,class privksP>
inline void CircuitBootstrappingFFT(TRGSWFFT<typename privksP::targetP> &trgswfft, const TLWE<typename bkP::domainP> &tlwe,
                             const CircuitKey<bkP,privksP> &ck))
{
    TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<bkP,privksP>(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * privksP::targetP::l; i++)
        for (int j = 0; j < 2; j++) TwistIFFT<privksP::targetP>(trgswfft[i][j], trgsw[i][j]);
}

void CircuitBootstrappingFFTlvl01(TRGSWFFT<lvl1param> &trgswfft, const TLWE<lvl0param> &tlwe,
                             const CircuitKey<lvl02param,lvl21param> &ck)
{
    CircuitBootstrappingFFT<lvl02param,lvl21param>(trgswfft,tlwe,ck);
}

void CircuitBootstrappingFFTlvl02(TRGSWFFT<lvl2param> &trgswfft, const TLWE<lvl0param> &tlwe,
                                const CircuitKey<lvl02param,lvl22param> &ck)
{
    CircuitBootstrappingFFT<lvl02param,lvl22param>(trgswfft,tlwe,ck);
}

void CircuitBootstrappingFFTInv(TRGSWFFTlvl1 &invtrgswfft, const TLWElvl0 &tlwe,
                                const CircuitKey &ck)
{
    TLWElvl0 invtlwe;
    // HomNot
    for (int i = 0; i <= DEF_n; i++) invtlwe[i] = -tlwe[i];
    TRGSWlvl1 trgsw;
    CircuitBootstrapping(trgsw, invtlwe, ck);
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFTlvl1(invtrgswfft[i][j], trgsw[i][j]);
}

void CircuitBootstrappingFFTwithInv(TRGSWFFTlvl1 &trgswfft,
                                    TRGSWFFTlvl1 &invtrgswfft,
                                    const TLWElvl0 &tlwe, const CircuitKey &ck)
{
    TRGSWlvl1 trgsw;
    array<uint32_t, DEF_l> h;
    for (int i = 0; i < DEF_l; i++) h[i] = 1U << (32 - (i + 1) * DEF_Bgbit);
    CircuitBootstrapping(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl1(trgswfft[i][j], trgsw[i][j]);
    for (int i = 0; i < DEF_l; i++) {
        for (int j = 0; j < DEF_N; j++) {
            trgsw[i][0][j] *= -1;
            trgsw[i][1][j] *= -1;
            trgsw[i + DEF_l][0][j] *= -1;
            trgsw[i + DEF_l][1][j] *= -1;
        }
        trgsw[i][0][0] += h[i];
        trgsw[i + DEF_l][1][0] += h[i];
    }
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFTlvl1(invtrgswfft[i][j], trgsw[i][j]);
}
}  // namespace TFHEpp