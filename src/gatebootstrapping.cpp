#include <cloudkey.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;
template <typename T = uint32_t, uint32_t N = DEF_N>
inline void PolynomialMulByXai(array<T, N> &res, const array<T, N> &poly,
                               const T a)
{
    if (a == 0)
        return;
    else if (a < N) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + N];
        for (int i = a; i < N; i++) res[i] = poly[i - a];
    }
    else {
        const T aa = a - N;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + N];
        for (int i = aa; i < N; i++) res[i] = -poly[i - aa];
    }
}

void PolynomialMulByXailvl1(Polynomiallvl1 &res, const Polynomiallvl1 &poly,
                            const uint32_t a)
{
    PolynomialMulByXai<uint32_t, DEF_N>(res, poly, a);
}

void PolynomialMulByXailvl2(Polynomiallvl2 &res, const Polynomiallvl2 &poly,
                            const uint64_t a)
{
    PolynomialMulByXai<uint64_t, DEF_nbar>(res, poly, a);
}

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void PolynomialMulByXaiMinusOne(array<T, N> &res,
                                       const array<T, N> &poly, const T a)
{
    if (a == 0)
        return;
    else if (a < N) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + N] - poly[i];
        for (int i = a; i < N; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const T aa = a - N;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + N] - poly[i];
        for (int i = aa; i < N; i++) res[i] = -poly[i - aa] - poly[i];
    }
}

void PolynomialMulByXaiMinusOnelvl1(Polynomiallvl1 &res,
                                    const Polynomiallvl1 &poly,
                                    const uint32_t a)
{
    PolynomialMulByXaiMinusOne<uint32_t, DEF_N>(res, poly, a);
}

inline void PolynomialMulByXaiMinusOnelvl2(Polynomiallvl2 &res,
                                           const Polynomiallvl2 &poly,
                                           const uint64_t a)
{
    PolynomialMulByXaiMinusOne<uint64_t, DEF_nbar>(res, poly, a);
}

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void RotatedTestVector(array<array<T, N>, 2> &testvector, uint32_t bara,
                              const T μ)
{
    testvector[0] = {};
    if (bara < N) {
        for (int i = 0; i < bara; i++) testvector[1][i] = -μ;
        for (int i = bara; i < N; i++) testvector[1][i] = μ;
    }
    else {
        const T baraa = bara - N;
        for (int i = 0; i < baraa; i++) testvector[1][i] = μ;
        for (int i = baraa; i < N; i++) testvector[1][i] = -μ;
    }
}

void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWElvl1 &acc, const TLWElvl0 &tlwe,
                                         const GateKey &gk)
{
    TRLWElvl1 temp;
    uint32_t bara = 2 * DEF_N - modSwitchFromTorus32<2 * DEF_N>(tlwe[DEF_n]);
    RotatedTestVector<uint32_t, DEF_N>(acc, bara, DEF_μ);
    for (int i = 0; i < DEF_n; i++) {
        bara = modSwitchFromTorus32<2 * DEF_N>(tlwe[i]);
        if (bara == 0) continue;
        PolynomialMulByXaiMinusOnelvl1(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOnelvl1(temp[1], acc[1], bara);
        trgswfftExternalProductlvl1(temp, temp, gk.bkfftlvl01[i]);
        for (int i = 0; i < DEF_N; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
}

void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res, const TLWElvl0 &tlwe,
                                        const GateKey &gk)
{
    TRLWElvl1 acc;
    GateBootstrappingTLWE2TRLWEFFTlvl01(acc, tlwe, gk);
    SampleExtractIndexlvl1(res, acc, 0);
}

inline void GateBootstrappingTLWE2TLWEFFTlvl02(TLWElvl2 &res,
                                               const TLWElvl0 &tlwe,
                                               const CircuitKey &ck,
                                               const uint64_t μs2)
{
    TRLWElvl2 acc;
    TRLWElvl2 temp;
    uint32_t bara =
        2 * DEF_nbar - modSwitchFromTorus64<2 * DEF_nbar>(tlwe[DEF_n]);
    RotatedTestVector<uint64_t, DEF_nbar>(acc, bara, μs2);
    for (int i = 0; i < DEF_n; i++) {
        bara = modSwitchFromTorus64<2 * DEF_nbar>(tlwe[i]);
        if (bara == 0) continue;
        PolynomialMulByXaiMinusOnelvl2(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOnelvl2(temp[1], acc[1], bara);
        trgswfftExternalProductlvl2(temp, temp, ck.bkfftlvl02[i]);
        for (int i = 0; i < DEF_nbar; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
    SampleExtractIndexlvl2(res, acc, 0);
    res[DEF_nbar] += μs2;
}

void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, const GateKey &gk)
{
    TLWElvl1 tlwelvl1;
    GateBootstrappingTLWE2TLWEFFTlvl01(tlwelvl1, tlwe, gk);
    IdentityKeySwitchlvl10(res, tlwelvl1, gk.ksk);
}

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