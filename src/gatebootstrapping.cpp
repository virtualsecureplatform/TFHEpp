#include <cloudkey.hpp>
#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>
#include <detwfa.hpp>

namespace TFHEpp {
using namespace std;

template <class P>
inline void GateBootstrappingTLWE2TRLWEFFT(
    TRLWE<typename P::targetP> &acc, const TLWE<typename P::domainP> &tlwe,
    const BootStrappingKeyFFT<P> &bkfft)
{
    TRLWE<typename P::targetP> temp;
    uint32_t bara = 2 * P::targetP::n - modSwitchFromTorus<typename P::targetP>(
                                            tlwe[P::domainP::n]);
    RotatedTestVector<typename P::targetP>(acc, bara, P::targetP::μ);
    for (int i = 0; i < P::domainP::n; i++) {
        bara = modSwitchFromTorus<typename P::targetP>(tlwe[i]);
        if (bara == 0) continue;
        // Do not use CMUX to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOne<typename P::targetP>(acc,bkfft[i],bara);
    }
}

template <class P>
inline void GateBootstrappingTLWE2TLWEFFT(TLWE<typename P::targetP> &res,
                                          const TLWE<typename P::domainP> &tlwe,
                                          const BootStrappingKeyFFT<P> &bkfft)
{
    TRLWE<typename P::targetP> acc;
    GateBootstrappingTLWE2TRLWEFFT<P>(acc, tlwe, bkfft);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

void GateBootstrappingTLWE2TRLWEFFTlvl01(
    TRLWE<lvl1param> &acc, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl01param> &bkfft)
{
    GateBootstrappingTLWE2TRLWEFFT<lvl01param>(acc, tlwe, bkfft);
}

void GateBootstrappingTLWE2TLWEFFTlvl01(
    TLWE<lvl1param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl01param> &bkfft)
{
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwe, bkfft);
}

void GateBootstrappingTLWE2TLWEFFTvariableMulvl02(
    TLWE<lvl2param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl02param> &bkfftlvl02, const lvl2param::T μs2)
{
    GateBootstrappingTLWE2TLWEFFTvariableMu<lvl02param>(res, tlwe, bkfftlvl02,
                                                        μs2);
}

void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const GateKey &gk)
{
    TLWE<lvl1param> tlwelvl1;
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(tlwelvl1, tlwe, gk.bkfftlvl01);
    IdentityKeySwitch<lvl10param>(res, tlwelvl1, gk.ksk);
}
}  // namespace TFHEpp