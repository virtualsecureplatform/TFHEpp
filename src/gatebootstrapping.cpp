#include<gatebootstrapping.hpp>
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
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a];
    }
    else {
        const typename P::T aa = a - P::n;
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

template<class P>
inline void GateBootstrappingTLWE2TRLWEFFT(TRLWE<typename P::targetP> &acc, const TLWE<typename P::domainP> &tlwe,
                                         const BootStrappingKeyFFT<P> &bkfft)
{
    TRLWE<typename P::targetP> temp;
    uint32_t bara = 2 * P::targetP::n - modSwitchFromTorus<typename P::targetP>(tlwe[P::domainP::n]);
    RotatedTestVector<typename P::targetP>(acc, bara, P::targetP::μ);
    for (int i = 0; i < P::domainP::n; i++) {
        bara = modSwitchFromTorus<typename P::targetP>(tlwe[i]);
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

template<class P>
inline void GateBootstrappingTLWE2TLWEFFT(TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
                                        const BootStrappingKeyFFT<P> &bkfft)
{
    TRLWE<typename P::targetP> acc;
    GateBootstrappingTLWE2TRLWEFFT<P>(acc, tlwe, bkfft);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWE<lvl1param> &acc, const TLWE<lvl0param> &tlwe,
                                         const BootStrappingKeyFFT<lvl01param> &bkfft)
{
    GateBootstrappingTLWE2TRLWEFFT<lvl01param>(acc,tlwe,bkfft);
}

void GateBootstrappingTLWE2TLWEFFTlvl01(TLWE<lvl1param> &res, const TLWE<lvl0param> &tlwe,
                                        const BootStrappingKeyFFT<lvl01param> &bkfft)
{
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(res,tlwe,bkfft);
}

void GateBootstrappingTLWE2TLWEFFTvariableMulvl02(
    TLWE<lvl2param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl02param> &bkfftlvl02, const lvl2param::T μs2)
{
    GateBootstrappingTLWE2TLWEFFTvariableMu<lvl02param>(res,tlwe,bkfftlvl02,μs2);
}

void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe, const GateKey &gk)
{
    TLWE<lvl1param> tlwelvl1;
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(tlwelvl1, tlwe, gk.bkfftlvl01);
    IdentityKeySwitch<lvl10param>(res, tlwelvl1, gk.ksk);
}
}  // namespace TFHEpp