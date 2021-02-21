#include <gatebootstrapping.hpp>
#include <mulfft.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {
using namespace std;

inline void CircuitBootstrapping(TRGSWlvl1 &trgsw, const TLWElvl0 &tlwe,
                                 const CircuitKey &ck)
{
    TLWElvl2 tlwelvl2;
    for (int i = 0; i < DEF_l; i++) {
        GateBootstrappingTLWE2TLWEFFTlvl02(
            tlwelvl2, tlwe, ck, 1UL << (64 - (i + 1) * DEF_Bgbit - 1));
        PrivKeySwitchlvl21(trgsw[i], tlwelvl2, 0, ck.privksk);
        PrivKeySwitchlvl21(trgsw[i + DEF_l], tlwelvl2, 1, ck.privksk);
    }
}

void CircuitBootstrappingFFT(TRGSWFFTlvl1 &trgswfft, const TLWElvl0 &tlwe,
                             const CircuitKey &ck)
{
    TRGSWlvl1 trgsw;
    CircuitBootstrapping(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl1(trgswfft[i][j], trgsw[i][j]);
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