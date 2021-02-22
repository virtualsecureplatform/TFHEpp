#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void PolynomialMulByXaiMinusOnelvl1(Polynomiallvl1 &res,
                                    const Polynomiallvl1 &poly,
                                    const uint32_t a);
void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWElvl1 &acc, const TLWElvl0 &tlwe,
                                         const GateKey &gk);
void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res, const TLWElvl0 &tlwe,
                                        const GateKey &gk);
void GateBootstrappingTLWE2TLWEFFTlvl02(
    TLWElvl2 &res, const TLWElvl0 &tlwe,
    const BootStrappingKeyFFTlvl02 &bkfftlvl02, const uint64_t μs2);
void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, const GateKey &gk);
}  // namespace TFHEpp