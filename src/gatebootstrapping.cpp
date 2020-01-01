#include <cloudkey.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

#include <iostream>

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

inline void PolynomialMulByXaiMinusOnelvl1(Polynomiallvl1 &res,
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

inline void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res,
                                               const TLWElvl0 &tlwe,
                                               const CloudKey &ck)
{
    TRLWElvl1 acc;
    TRLWElvl1 temp;
    uint32_t bara = 2 * DEF_N - modSwitchFromTorus32<2 * DEF_N>(tlwe[DEF_n]);
    RotatedTestVector<uint32_t, DEF_N>(acc, bara, DEF_μ);
    for (int i = 0; i < DEF_n; i++) {
        bara = modSwitchFromTorus32<2 * DEF_N>(tlwe[i]);
        if (bara == 0) continue;
        PolynomialMulByXaiMinusOnelvl1(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOnelvl1(temp[1], acc[1], bara);
        trgswfftExternalProductlvl1(temp, ck.bkfftlvl01[i]);
        for (int i = 0; i < DEF_N; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
    SampleExtractIndexlvl1(res, acc, 0);
}

void GateBootstrappingTLWE2TLWEFFTlvl02(TLWElvl2 &res, const TLWElvl0 &tlwe,
                                        const CloudKey &ck, const uint64_t μs2)
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
        trgswfftExternalProductlvl2(temp, ck.bkfftlvl02[i]);
        for (int i = 0; i < DEF_nbar; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
    SampleExtractIndexlvl2(res, acc, 0);
    res[DEF_nbar] += μs2;
}

void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, const CloudKey &ck)
{
    TLWElvl1 tlwelvl1;
    GateBootstrappingTLWE2TLWEFFTlvl01(tlwelvl1, tlwe, ck);
    IdentityKeySwitchlvl10(res, tlwelvl1, ck.ksk);
}

void CircuitBootstrapping(TRGSWlvl1 &trgsw, const TLWElvl0 &tlwe,
                          const CloudKey &ck)
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
                             const CloudKey &ck)
{
    TRGSWlvl1 trgsw;
    CircuitBootstrapping(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl1(trgswfft[i][j], trgsw[i][j]);
}
}  // namespace TFHEpp