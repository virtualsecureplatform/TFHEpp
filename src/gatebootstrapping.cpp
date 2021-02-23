#include <cloudkey.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;
template <class P>
inline void PolynomialMulByXai(array<typename P::T, P::n> &res, const array<typename P::T, P::n> &poly,
                               const typename P::T a)
{
    if (a == 0)
        return;
    else if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + N];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a];
    }
    else {
        const P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa];
    }
}

void PolynomialMulByXailvl1(Polynomial<lvl1param> &res, const Polynomial<lvl1param> &poly,
                            const lvl1param::T a)
{
    PolynomialMulByXai<lvl1param>(res, poly, a);
}

void PolynomialMulByXailvl2(Polynomial<lvl2param> &res, const Polynomial<lvl2param> &poly,
                            const lvl2param::T a)
{
    PolynomialMulByXai<lvl2param>(res, poly, a);
}

template <class P>
inline void PolynomialMulByXaiMinusOne(array<typename P::T, P::n> &res,
                                       const array<typename P::T, P::n> &poly, const typename P::T a)
{
    if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n] - poly[i];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n] - poly[i];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa] - poly[i];
    }
}

void PolynomialMulByXaiMinusOnelvl1(Polynomial<lvl1param> &res,
                                    const Polynomial<lvl1param> &poly,
                                    const lvl1param::T a)
{
    PolynomialMulByXaiMinusOne<lvl1param>(res, poly, a);
}

inline void PolynomialMulByXaiMinusOnelvl2(Polynomial<lvl2param> &res,
                                           const Polynomial<lvl2param> &poly,
                                           const lvl2param::T a)
{
    PolynomialMulByXaiMinusOne<lvl2param>(res, poly, a);
}

template <class P>
inline void RotatedTestVector(array<array<typename P::T, P::n>, 2> &testvector, uint32_t bara,
                              const typename P::T μ)
{
    testvector[0] = {};
    if (bara < P::n) {
        for (int i = 0; i < bara; i++) testvector[1][i] = -μ;
        for (int i = bara; i < P::n; i++) testvector[1][i] = μ;
    }
    else {
        const P::T baraa = bara - P::n;
        for (int i = 0; i < baraa; i++) testvector[1][i] = μ;
        for (int i = baraa; i < P::n; i++) testvector[1][i] = -μ;
    }
}

template<class P>
inline void GateBootstrappingTLWE2TRLWEFFT(TRLWE<typename P::targetP> &acc, const TLWE<typename P::domainP> &tlwe,
                                         const BootStrappingKeyFFT<P> &bkfft)
{
    TRLWE<typename P::targetP> temp;
    uint32_t bara = 2 * P::targetP::n - modSwitchFromTorus<P::targetP>(tlwe[P::domainP::n]);
    RotatedTestVector<P::targetP>(acc, bara, P::targetP::μ);
    for (int i = 0; i < P::domainP::n; i++) {
        bara = modSwitchFromTorus<P::targetP>(tlwe[i]);
        if (bara == 0) continue;
        // Do not use CMUX to avoid unnecessary copy.
        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[1], acc[1], bara);
        trgswfftExternalProduct<typename P::targetP>(temp, temp, bkfft[i]);
        for (int i = 0; i < P::targetP::n; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
}

void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWE<lvl1param> &acc, const TLWE<lvl0param> &tlwe,
                                         const BootStrappingKeyFFT<lvl01param> &bkfft)
{
    GateBootstrappingTLWE2TRLWEFFT<lvl01param>(acc,tlwe,bkfft);
}

void GateBootstrappingTLWE2TLWEFFTlvl01(TLWE<lvl1param> &res, const TLWE<lvl0param> &tlwe,
                                        const BootStrappingKeyFFT<lvl01param> &bkfft)
{
    TRLWE<lvl1param> acc;
    GateBootstrappingTLWE2TRLWEFFTlvl01(acc, tlwe, bkfft);
    SampleExtractIndexlvl1(res, acc, 0);
}

template<class P>
void GateBootstrappingTLWE2TLWEFFTvaribaleMu(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootStrappingKeyFFT<P> &bkfft, const typename P::targetP::T μs2)
{
    TRLWE<typename P::targetP> acc,temp;
    uint32_t bara =
        2 * P::targetP::nbar - modSwitchFromTorus<typename P::targetP>(tlwe[P::domainP::n]);
    RotatedTestVector<P::targetP::T, P::targetP::n>(acc, bara, μs2);
    for (int i = 0; i < P::domainP::n; i++) {
        bara = modSwitchFromTorus<typename P::targetP>(tlwe[i]);
        if (bara == 0) continue;
        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[1], acc[1], bara);
        trgswfftExternalProduct<typename P::targetP>(temp, temp, bkfftlvl02[i]);
        for (int i = 0; i < DEF_nbar; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }
    SampleExtractIndexlvl2(res, acc, 0);
    res[DEF_nbar] += μs2;
}

void GateBootstrappingTLWE2TLWEFFTlvl02(
    TLWE<lvl2param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl02param> &bkfftlvl02, const lvl2param::T μs2)
{
    TRLWE<lvl2param> acc;
    TRLWElvl2 temp;
    uint32_t bara =
        2 * DEF_nbar - modSwitchFromTorus64<2 * DEF_nbar>(tlwe[DEF_n]);
    RotatedTestVector<uint64_t, DEF_nbar>(acc, bara, μs2);
    for (int i = 0; i < DEF_n; i++) {
        bara = modSwitchFromTorus64<2 * DEF_nbar>(tlwe[i]);
        if (bara == 0) continue;
        PolynomialMulByXaiMinusOnelvl2(temp[0], acc[0], bara);
        PolynomialMulByXaiMinusOnelvl2(temp[1], acc[1], bara);
        trgswfftExternalProductlvl2(temp, temp, bkfftlvl02[i]);
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
}  // namespace TFHEpp