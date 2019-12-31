#include <cloudkey.hpp>
#include <keyswitch.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
using namespace std;
template <typename T = uint32_t, uint32_t N = DEF_N>
inline void PolynomialMulByXai(array<T, N> &res, array<T, N> &poly, const T a)
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
void PolynomialMulByXailvl1(Polynomiallvl1 &res, Polynomiallvl1 &poly,
                            const uint32_t a)
{
    PolynomialMulByXai<uint32_t, DEF_N>(res, poly, a);
}

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void PolynomialMulByXaiMinusOne(array<T, N> &res, array<T, N> &poly,
                                       const T a)
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
                                           Polynomiallvl1 &poly,
                                           const uint32_t a)
{
    PolynomialMulByXaiMinusOne<uint32_t, DEF_N>(res, poly, a);
}

template <typename T = uint32_t, uint32_t N = DEF_N, T MU = DEF_MU>
inline void RotatedTestVector(array<array<T, N>, 2> &testvector, uint32_t bara)
{
    testvector[0] = {};
    if (bara < N) {
        for (int i = 0; i < bara; i++) testvector[1][i] = -MU;
        for (int i = bara; i < N; i++) testvector[1][i] = MU;
    }
    else {
        const T baraa = bara - N;
        for (int i = 0; i < baraa; i++) testvector[1][i] = MU;
        for (int i = baraa; i < N; i++) testvector[1][i] = -MU;
    }
}

inline void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res,
                                               const TLWElvl0 &tlwe,
                                               CloudKey &ck)
{
    TRLWElvl1 acc;
    TRLWElvl1 temp;
    uint32_t bara = 2 * DEF_N - (tlwe[DEF_n] >> (32 - (DEF_Nbit + 1)));
    RotatedTestVector<uint32_t, DEF_N, DEF_MU>(acc, bara);
    for (int i = 0; i < DEF_n; i++) {
        bara = (tlwe[i] >> (32 - (DEF_Nbit + 1)));
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

void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, CloudKey &ck)
{
    TLWElvl1 tlwelvl1;
    GateBootstrappingTLWE2TLWEFFTlvl01(tlwelvl1, tlwe, ck);
    IdentityKeySwitchlvl10(res, tlwelvl1, ck.ksk);
}
}  // namespace TFHEpp